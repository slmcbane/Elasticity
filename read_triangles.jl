using WriteVTK

read_triangles(; input::AbstractString, output::AbstractString) = open(input) do f
    num_nodes = read(f, Int)
    num_elements = read(f, Int)
    coords = Matrix{Float64}(undef, 2, num_nodes)
    read!(f, coords)
    vertices = Matrix{Int64}(undef, 3, num_elements)
    read!(f, vertices)
    vertices .+= 1

    cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, vertices[:, j])
             for j in 1:num_elements]
    
    vtk_grid(output, coords, cells)
end
