function dnf=funct_finiteDer_2D(x,coeff)
    % This function simulates the diffusion of the activator and/or inhibitor (when applicable)
    dnf=imfilter(x,coeff);
end