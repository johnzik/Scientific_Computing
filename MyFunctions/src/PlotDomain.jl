function PlotDomain(all_coords::Matrix{Float64})
    fig_mesh = GLMakie.Figure();

    # Plot domain with subdomain 1 overlaps
    ax_mesh = GLMakie.Axis(fig_mesh[1, 1], aspect = AxisAspect(1), xlabel = "x", ylabel = "y", 
                           title = "Physical Mesh Nodes of L-Shape Domain");

    GLMakie.scatter!(ax_mesh, all_coords[:,1], all_coords[:,2],
                    color = :blue, markersize = 5,
                    strokecolor = :black, strokewidth = 0.5, # Stroke is for better visibility
                    );
    
    display(GLMakie.Screen(), fig_mesh);
end