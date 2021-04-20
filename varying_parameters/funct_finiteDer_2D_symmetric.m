function dnf=funct_finiteDer_2D_symmetric(x,coeff)
    % This function simulates the diffusion of the activator when varying parameters 
    % using a domain which generates planar waves, thus we use reflecting boundary conditions
    dnf=imfilter(x,coeff,'symmetric');
end