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