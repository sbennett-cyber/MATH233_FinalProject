% Create the parameter files

NumChannels = input("How many edges are in the network? Options are 1, 2, or 3?")
NumFields   = input("How many fields are in the domain? Options are 1, 2, or 3?")

Channels = input( "Provide the coordinates and ID number for each edge of the network as a N x 3 matrix with the last column listing the ID number starting with 1.")

switch NumChannels
    case 1
    case 2
    case 3
end 

