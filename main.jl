


using Random, Distributions

include("functions.jl")



#============================================================

8/11/2018

Tumour vector reformulated as a matrix, where a cell with 
coordinates (x,y) will occupy the [x,y] element of matrix.


============================================================#

#=
# Define a cell 
mutable struct cell
	
	dvr::Int32
	res::Int32
	pgr::Int32
	#posX::Int32
	#posY::Int32

end


# Function which finds index of cell in tumour positioned at given coordinates
function findposition(tumour , x::Int , y::Int)

	# Construct two arrays with indices of cells in tumour at given x and y values
	xhits = findall(in([x]) , getfield.(tumour, :posX))
	yhits = findall(in([y]) , getfield.(tumour, :posY))

	# Index of cell we want is that which is common to both arrays 
	index = intersect(xhits , yhits)
	
	return index 			# Returns array with indices of common elements in xhits and yhits (so will be a 1d array with a single entry)

end



# Define cell division
function divide(tumour , cell_index_x::Int , cell_index_y::Int , x::Int , y::Int , tumour_size::Int)

	queue = 0

	# Count how many cells need to be pushed in the specified direction (quantified by queue variable)
	while true

		#filled = false
		if ( tumour[cell_index_x + x*(queue+1) , cell_index_y + y*(queue+1)].dvr == -1 ) # if not occupied
			break
		end

		queue += 1

	end 


	# Shove glands outwards
	for j in 1:queue

		tumour[cell_index_x + x*(queue - j + 2) , cell_index_y + y*(queue - j + 2)] = tumour[cell_index_x + x*(queue - j + 1) , cell_index_y + y*(queue - j + 1)]

	end


	# Create daughter cell
	tumour[cell_index_x + x , cell_index_y + y] = cell(tumour[cell_index_x , cell_index_y].dvr , tumour[cell_index_x, cell_index_y].res , tumour[cell_index_x , cell_index_y].pgr)

	tumour_size += 1

	return tumour_size

end


# Round and parse a given floating point to integer
function float_to_int(x)

	minval = min(x - floor(x) , abs(x - ceil(x)))
	y = floor(x) + (findfirst(isequal(minval) , [x - floor(x) , abs(x - ceil(x))]) - 1)
	y_int = trunc(Int , y)

	return y_int
    
end
=#


#================== Create output files ===================#



# Check for existing ./DATA directory, and create if it does not exist 
if !isdir("./DATA") mkdir("./DATA") end


# Open data file
filename = "DATA/"
open("params.dat") do file

	for ln in eachline(file)
		global filename *= ln
		global filename *= "_"
	end

end


# Open file with number of cells versus time 
filename = chop(filename)
filename *= "_Nversust.dat"
NversusT_file = open(filename , "w")





#================= Process steering file ====================#



# Read simulation paramaters from steering file ("params.dat")
println("----------------------")

line = 1
params = []
open("params.dat") do file
	for ln in eachline(file)
		println(ln)
		index = findfirst(isequal('=') , ln)
		if (!occursin(".", ln)) push!(params , parse(Int32 , ln[index+1:end])) # If argument is intended to be an integer
		elseif (occursin(".", ln)) push!(params , parse(Float64 , ln[index+1:end])) # If argument is intended to be a floating point 
		end
		global line += 1
	end
end

println("----------------------")



# Seed the random number generator
Random.seed!(params[2])



#================== Initialise tumour ====================#

# Estimate radius of resulting tumour using fitted parameters from previous simulations (slightly over-estimate)
#area = 10.0*exp(0.68572*params[1])
radius = (params[1]/pi)^0.5
radius = float_to_int(1.25*radius)


# Define tumour as an array of cells ("normal" tissue modelled as cell with all GAs = -1, ignored by algorithm)
tumour = fill( cell(-1,-1,-1) , 2*radius , 2*radius )

println("Lattice size: $(2*radius)")


println("")
println("Initialising tumour... Done")
println("")


# Seed first tumour cell/gland at (x,y) = (0,0)
tumour[radius , radius] = cell(0,0,0)
global tumour_size = 1


# Set birth/death rate ratio
const bdratio = 0.5

#= 
	Define mutation probailities for driver (d), resistant (r) and passenger 
   	mutations (latter defined implicitly through total mutation probability (t))
=#

# Define Poisson distributions
poisson_t = Poisson(params[4] - params[5] - params[6])
poisson_d = Poisson(params[5])
poisson_r = Poisson(params[6])


# Reset time variable
t = 0.0


const r_death = bdratio * log(2.0)




#================== Simulate tumour growth ==================#

#inserted = false
iter = 0
max_birth = 0.0
x = 0
y = 0
x_tot = 0
y_tot = 0
while (true)

	global iter += 1

	# Compute estimate of tumour "radius"
	rad = ceil((tumour_size/pi)^0.5)
	rad = float_to_int(rad)

 
	# Randomly select one cell to divide
	while(true)
		global cell_index_x = rand(radius-rad:radius+rad)
		global cell_index_y = rand(radius-rad:radius+rad)
		#println("$cell_index_x $cell_index_x")
		if !(tumour[cell_index_x , cell_index_y].dvr == -1) break end
	end


	# Compute birth and death rate of cell (params[3] is the selective advantage of a single driver mutation)
	r_birth = log(2.0) * ((1.0 + params[3])^tumour[cell_index_x , cell_index_y].dvr)
	#r_death = bdratio * r_birth 		# Model 1
	#r_death = bdratio * log(2.0)		# Model 2


	# Update maximal birth and death rate of all cells 
	if (r_birth > max_birth) max_birth = r_birth end


	# Randomly select the direction in which to divide
	while(true)
           global x = rand(-1:1)
           global y = rand(-1:1)
           if !((x == 0) && (y == 0)) break end
    end


	# Cell divides with proability r_birth/max_birth
	if (rand() < (r_birth/max_birth))

		# Simulate cell division
		tumour_size = divide(tumour , cell_index_x , cell_index_y , x , y , tumour_size)

		# Add new GAs to daughter cells
		#initial_dvr = tumour[cell_index_x , cell_index_y].dvr

		tumour[cell_index_x , cell_index_y].dvr += rand(poisson_d)
		tumour[cell_index_x , cell_index_y].res += rand(poisson_r)
		tumour[cell_index_x , cell_index_y].pgr += rand(poisson_t)

		#final_dvr = tumour[cell_index_x , cell_index_y].dvr
		#if (final_dvr-initial_dvr != 0) println("New driver GA at (x,y) = ($cell_index_x,$cell_index_y)") end


		tumour[cell_index_x + x , cell_index_y + y].dvr += rand(poisson_d)
		tumour[cell_index_x + x , cell_index_y + y].res += rand(poisson_r)
		tumour[cell_index_x + x , cell_index_y + y].pgr += rand(poisson_t)


		#if (tumour[cell_index].dvr == 1 || tumour[tumour_size].dvr == 1) println("New driver mutation after n=$iter steps!") end


	elseif (rand() < (r_death/max_birth))
 	
		# Delete cell from tumour
		tumour[cell_index_x , cell_index_y] = cell(-1,-1,-1)

		# Size of tumour is reduced by 1
		global tumour_size -= 1

	end

	# Progress time variable
	global t += 1.0/(max_birth * tumour_size)

	# Write total number of cells after regular number of iterations
	if (iter%200 == 0)

		#max_dvr = log(max_birth)/log(log(2.0) * ((1.0 + params[3])))
		#max_dvr = float_to_int(max_dvr)			# Find maximal number of drivers in an individual cell in tumour

		println(NversusT_file , "$t $tumour_size")
		println("Iter=$iter, N=$tumour_size")

		# Update max_birth variable if needs be
		max_dvr = findmax(getfield.(tumour , :dvr))[1]			# getfield.(tumour, :dvr) returns array dvr value of all cells in tumour
		global max_birth = log(2.0) * ((1.0 + params[3])^(max_dvr))
		#println("Updated max_birth variable to $max_birth \n")

	end


	# Update global tumour size variable 
	#global tumour_size = current_size


	if (tumour_size > params[1]) break end


end

#println("<x> = $x_tot")
#println("<y> = $y_tot")




#================== Write data to files ====================#


# Open data file
filename = "DATA/"
open("params.dat") do file

	for ln in eachline(file)
		global filename *= ln
		global filename *= "_"
	end

end

# Open main data file containing tumour
filename = chop(filename)
filename *= ".dat"
tumourfile = open(filename , "w")


# Write data
for i in 1:2*radius
	for j in 1:2*radius
		if (tumour[i,j].dvr != -1)
			println(tumourfile , "$i  $j $(tumour[i,j].dvr) $(tumour[i,j].res) $(tumour[i,j].pgr)")
		end
	end
end





