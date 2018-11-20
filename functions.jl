



using Random, Distributions



# Define a cell 
mutable struct cell
	
	dvr::Int32
	res::Int32
	pgr::Int32

end


# Read input partameters from steering file
line = 1
params = []
open("params.dat") do file
	for ln in eachline(file)
		index = findfirst(isequal('=') , ln)
		if (!occursin(".", ln)) push!(params , parse(Int32 , ln[index+1:end])) # If argument is intended to be an integer
		elseif (occursin(".", ln)) push!(params , parse(Float64 , ln[index+1:end])) # If argument is intended to be a floating point 
		end
		global line += 1
	end
end

# Seed the random number generator
Random.seed!(params[2])

# Define Poisson distributions
poisson_t = Poisson(params[4] - params[5] - params[6])
poisson_d = Poisson(params[5])
poisson_r = Poisson(params[6])



# Function which finds index of cell in tumour positioned at given coordinates
function findposition(tumour , x::Int , y::Int)

	# Construct two arrays with indices of cells in tumour at given x and y values
	xhits = findall(in([x]) , getfield.(tumour, :posX))
	yhits = findall(in([y]) , getfield.(tumour, :posY))

	# Index of cell we want is that which is common to both arrays 
	index = intersect(xhits , yhits)
	
	return index 			# Returns array with indices of common elements in xhits and yhits (so will be a 1d array with a single entry)

end



# Define cell division in volumetric growth model
function MODEL1_divide(tumour , cell_index_x::Int , cell_index_y::Int, tumour_size::Int)

	# Randomly select the direction in which to divide
	while(true)
           global x = rand(-1:1)
           global y = rand(-1:1)
           if !((x == 0) && (y == 0)) break end
    end

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

	# Add new GAs to daughter cells
	tumour[cell_index_x , cell_index_y].dvr += rand(poisson_d)
	tumour[cell_index_x , cell_index_y].res += rand(poisson_r)
	tumour[cell_index_x , cell_index_y].pgr += rand(poisson_t)

	tumour[cell_index_x + x , cell_index_y + y].dvr += rand(poisson_d)
	tumour[cell_index_x + x , cell_index_y + y].res += rand(poisson_r)
	tumour[cell_index_x + x , cell_index_y + y].pgr += rand(poisson_t)

	return tumour_size

end



# Define cell division in volumetric growth model, but where cells completely surrounded by other tumour cells cannot replicate (p_divide NOT proportional to number of empty neighbours)
function MODEL2_divide(tumour, x_nn , y_nn , cell_index_x::Int , cell_index_y::Int , tumour_size::Int , empty_neighbours::Int)


	# If any empty neighbours, divide into one with flat probability
	if (empty_neighbours != 0)

		chosen_direction = rand(1:empty_neighbours)

		# Create daughter cell
		tumour[cell_index_x + x_nn[chosen_direction] , cell_index_y + y_nn[chosen_direction]] = cell(tumour[cell_index_x , cell_index_y].dvr , tumour[cell_index_x, cell_index_y].res , tumour[cell_index_x , cell_index_y].pgr)

		tumour_size += 1

		# Add new GAs to daughter cells
		tumour[cell_index_x , cell_index_y].dvr += rand(poisson_d)
		tumour[cell_index_x , cell_index_y].res += rand(poisson_r)
		tumour[cell_index_x , cell_index_y].pgr += rand(poisson_t)

		tumour[cell_index_x + x_nn[chosen_direction] , cell_index_y + y_nn[chosen_direction]].dvr += rand(poisson_d)
		tumour[cell_index_x + x_nn[chosen_direction] , cell_index_y + y_nn[chosen_direction]].res += rand(poisson_r)
		tumour[cell_index_x + x_nn[chosen_direction] , cell_index_y + y_nn[chosen_direction]].pgr += rand(poisson_t)
	end

	return tumour_size


end


# Define cell division in volumetric growth model, but where cells completely surrounded by other tumour cells cannot replicate (p_divide NOT proportional to number of empty neighbours)
function MODEL3_divide(tumour, cell_index_x::Int , cell_index_y::Int , tumour_size::Int)


	# Randomly select the direction in which to divide
	while(true)
           global x = rand(-1:1)
           global y = rand(-1:1)
           if !((x == 0) && (y == 0)) break end
    end

	if (tumour[cell_index_x + x , cell_index_y + y].dvr == -1)		# Check if neighbour is empty

		# Create daughter cell
		tumour[cell_index_x + x , cell_index_y + y] = cell(tumour[cell_index_x , cell_index_y].dvr , tumour[cell_index_x, cell_index_y].res , tumour[cell_index_x , cell_index_y].pgr)

		tumour_size += 1

		# Add new GAs to daughter cells
		tumour[cell_index_x , cell_index_y].dvr += rand(poisson_d)
		tumour[cell_index_x , cell_index_y].res += rand(poisson_r)
		tumour[cell_index_x , cell_index_y].pgr += rand(poisson_t)

		tumour[cell_index_x + x , cell_index_y + y].dvr += rand(poisson_d)
		tumour[cell_index_x + x , cell_index_y + y].res += rand(poisson_r)
		tumour[cell_index_x + x , cell_index_y + y].pgr += rand(poisson_t)

	end


	return tumour_size


end



# Round and parse a given floating point to integer
function float_to_int(x)

	minval = min(x - floor(x) , abs(x - ceil(x)))
	y = floor(x) + (findfirst(isequal(minval) , [x - floor(x) , abs(x - ceil(x))]) - 1)
	y_int = trunc(Int , y)

	return y_int
    
end








