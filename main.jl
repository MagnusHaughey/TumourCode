


using Random, Distributions

include("functions.jl")


# Declare global variables
const NEIGHBOURHOOD = 8
const bdratio = 0.5 	# birth/death-rate ratio



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


# Open file for number of cells versus time 
NversusT_filename = chop(filename)
NversusT_filename *= "_Nversust.dat"
NversusT_file = open(NversusT_filename , "w")

tumour_filename = chop(filename)
tumour_filename *= ".dat"



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


# Determine growth model
MODEL = fill(false , 3)		# List of booleans correpsonding to all available models
if ( (params[7] < 1) || (params[7] > length(MODEL)) )
	println("""Error with model variable in "params.dat" file. Exiting...""")
	exit(0)
end
MODEL[params[7]] = true



#================== Initialise tumour ====================#

# Estimate radius of resulting tumour using fitted parameters from previous simulations (slightly over-estimate)
radius = (params[1]/pi)^0.5
radius = float_to_int(1.25*radius)


# Define tumour as an array of cells ("normal" tissue modelled as cell with all GAs = -1, ignored by algorithm)
tumour = fill( cell(-1,-1,-1) , (2*radius , 2*radius) )

println("Lattice size: $(2*radius)")


println("")
println("Initialising tumour... Done")
println("")


# Seed first tumour cell/gland at (x,y) = (0,0)
tumour[radius , radius] = cell(0,0,0)
global tumour_size = 1


#= 
	Define mutation probailities for driver (d), resistant (r) and passenger 
   	mutations (latter defined implicitly through total mutation probability (t))
=#

# Define Poisson distributions
#poisson_t = Poisson(params[4] - params[5] - params[6])
#poisson_d = Poisson(params[5])
#poisson_r = Poisson(params[6])


# Reset time variable
t = 0.0


# Define constant death rate for all cells
const r_death = bdratio * log(2.0)		# Death model 1


# Define arrays which will contain relative coordinates of empty neighbours for a chosen cell (for model 2)
x_nn = fill(0 , NEIGHBOURHOOD)
y_nn = fill(0 , NEIGHBOURHOOD)



#================== Simulate tumour growth ==================#

#inserted = false
iter = 0
max_birth = 0.0
x = 0
y = 0

while (true)

	global iter += 1

	# Compute estimate of tumour "radius"
	#rad = ceil((tumour_size/pi)^0.5)
	#rad = float_to_int(1.1rad)			# Model 2 requires a slight over-estimation of radius to enable surface cells to be picked by algorithm
	
	#lower_bound = findmax([ 0 , radius-rad ])[1]

 
	# Randomly select one cell to divide
	while(true)
		#global cell_index_x = rand(lower_bound:radius+rad)
		#global cell_index_y = rand(lower_bound:radius+rad)

		global cell_index_x = rand(1:(2*radius))
		global cell_index_y = rand(1:(2*radius))
		if !(tumour[cell_index_x , cell_index_y].dvr == -1) break end
	end


	# Compute birth and death rate of cell (params[3] is the selective advantage of a single driver mutation)
	r_birth = log(2.0) * ((1.0 + params[3])^tumour[cell_index_x , cell_index_y].dvr)		# Birth model 1
	#r_birth = log(2.0) * (1.0 + ((tumour[cell_index_x , cell_index_y].dvr)*params[3]) )		# Birth model 2


	#r_death = bdratio * r_birth 		# Death model 2


	# Update maximal birth and death rate of all cells 
	if (r_birth > max_birth) max_birth = r_birth end



	# Cell divides with proability r_birth/max_birth
	if (rand() < (r_birth/max_birth))

		# Simulate cell division
		if MODEL[1]

			tumour_size = MODEL1_divide(tumour , cell_index_x , cell_index_y , tumour_size)

		elseif MODEL[2]

			empty_neighbours = 0

			# Check for any neighbouring empty lattice sites
			for i in -1:1
				for j in -1:1

					if (tumour[cell_index_x + i , cell_index_y + j].dvr == -1) 	# if not occupied
						empty_neighbours += 1
						x_nn[empty_neighbours] = i 	# Store coordinates of empty neighbour
						y_nn[empty_neighbours] = j
					end
				end
			end

			if (empty_neighbours != 0)
				tumour_size = MODEL2_divide(tumour, x_nn , y_nn , cell_index_x , cell_index_y , tumour_size, empty_neighbours)
			end

		elseif MODEL[3]

			tumour_size = MODEL3_divide(tumour , cell_index_x , cell_index_y , tumour_size)

		end


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

		println(NversusT_file , "$t $tumour_size")
		println("Iter=$iter, N=$tumour_size")

		# Update max_birth variable if needs be
		max_dvr = findmax(getfield.(tumour , :dvr))[1]			# getfield.(tumour, :dvr) returns array dvr value of all cells in tumour
		global max_birth = log(2.0) * ((1.0 + params[3])^(max_dvr))

	end


	if (tumour_size > params[1]) break end


end



#================== Write data to files ====================#



# Open main data file containing tumour
tumourfile = open(tumour_filename , "w")


# Write data
for i in 1:2*radius
	for j in 1:2*radius
		if (tumour[i,j].dvr != -1)
			println(tumourfile , "$i  $j $(tumour[i,j].dvr) $(tumour[i,j].res) $(tumour[i,j].pgr)")
		end
	end
end





