function PlotDomainOverlaps(all_coords::Matrix{Float64}, idx::Vector{Vector{Int64}})

    fig_mesh = GLMakie.Figure();

    # Plot domain with subdomain 1 overlaps
    ax_mesh = GLMakie.Axis(fig_mesh[2, 1:2], aspect = AxisAspect(1), xlabel = "x", ylabel = "y", 
                           title = "Physical Mesh Nodes of L-Shape Domain with Subdomain 1 Overlaps");

    GLMakie.scatter!(ax_mesh, all_coords[:,1], all_coords[:,2],
                    color = :blue, markersize = 5,
                    strokecolor = :black, strokewidth = 0.5, # Stroke is for better visibility
                    );

    GLMakie.scatter!(ax_mesh, all_coords[idx[1],1], all_coords[idx[1],2],
                    color = :red, markersize = 4,
                    strokecolor = :black, strokewidth = 0.5, # Stroke is for better visibility
                    );

    # Plot domain with subdomain 2 overlaps
    ax_mesh2 = GLMakie.Axis(fig_mesh[1,1], aspect = AxisAspect(1), xlabel = "x", ylabel = "y", 
                            title = "Physical Mesh Nodes of L-Shape Domain with Subdomain 2 Overlaps");

    GLMakie.scatter!(ax_mesh2, all_coords[:,1], all_coords[:,2],
                    color = :blue, markersize = 5,
                    strokecolor = :black, strokewidth = 0.5, # Stroke is for better visibility
                    );

    GLMakie.scatter!(ax_mesh2, all_coords[idx[2],1], all_coords[idx[2],2],
                    color = :red, markersize = 4,
                    strokecolor = :black, strokewidth = 0.5, # Stroke is for better visibility
                    );

    # Plot domain with subdomain 3 overlaps
    ax_mesh3 = GLMakie.Axis(fig_mesh[1, 2], aspect = AxisAspect(1), xlabel = "x", ylabel = "y", 
                            title = "Physical Mesh Nodes of L-Shape Domain with Subdomain 3 Overlaps");

    GLMakie.scatter!(ax_mesh3, all_coords[:,1], all_coords[:,2],
                    color = :blue, markersize = 5,
                    strokecolor = :black, strokewidth = 0.5, # Stroke is for better visibility
                    );

    GLMakie.scatter!(ax_mesh3, all_coords[idx[3],1], all_coords[idx[3],2],
                    color = :red, markersize = 4,
                    strokecolor = :black, strokewidth = 0.5, # Stroke is for better visibility
                    );

    display(GLMakie.Screen(), fig_mesh);
end