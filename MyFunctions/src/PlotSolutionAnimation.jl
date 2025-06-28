function PlotSolutionAnimation(GO::Vector{Matrix{Int64}}, u, k::Float64)
    GLMakie.activate!();
    rowNodes = length(GO[1][1,:]);
    newStep = 1/(rowNodes - 1); # Get step from the GO matrix

    c = 3*10^8; # Speed of the wave
    omega = k*c; # Angular Velocity
    T = 2π / omega; # Period
    frames = 60;
    ts = range(0, T, length=frames);  # Time frames

    # --- Mesh for plotting ---
    X1m, Y1m = ndgrid(-1:newStep:0, -1:newStep:0);
    X2m, Y2m = ndgrid(0:newStep:1, -1:newStep:0);
    X3m, Y3m = ndgrid(-1:newStep:0, 0:newStep:1);

    u1 = u[GO[1]]; U1 = reshape(u1, size(GO[1]));
    u2 = u[GO[2]]; U2 = reshape(u2, size(GO[2]));
    u3 = u[GO[3]]; U3 = reshape(u3, size(GO[3]));

    # --- Create Observable for each subdomain ---
    U1_t = Observable(real.(U1));
    U2_t = Observable(real.(U2));
    U3_t = Observable(real.(U3));

    # --- Create Figure ---
    fig2 = GLMakie.Figure();
    screen2 = display(GLMakie.Screen(), fig2);
    ax = GLMakie.Axis3(fig2[1,1], title = "Wave Propagation", xlabel="x", ylabel="y", zlabel="u(x,y,t)");

    # --- Surface Plots ------------------------------------
    GLMakie.surface!(ax, X1m, Y1m, U1_t, colormap = :viridis);
    GLMakie.surface!(ax, X2m, Y2m, U2_t, colormap = :viridis);
    GLMakie.surface!(ax, X3m, Y3m, U3_t, colormap = :viridis);

    # --- Animation ----------------------------------------
    for t in ts
        phase = exp(1im * omega * t);
        U1_t[] = real.(U1 * phase);
        U2_t[] = real.(U2 * phase);
        U3_t[] = real.(U3 * phase);
        sleep(0.05)
    end
end