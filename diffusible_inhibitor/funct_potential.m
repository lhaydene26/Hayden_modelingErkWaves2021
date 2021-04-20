function [VE,VA,VI]=funct_potential(E,A,I,params)
    % This function calculates the reaction terms in the reaction-diffusion models

    alpha_1 = params.alpha_1;
    beta_1 = params.beta_1;
    gamma_1 = params.gamma_1;
    alpha_2 = params.alpha_2;
    alpha_3 = params.alpha_3;
    gamma_2 = params.gamma_2;
    alpha_4 = params.alpha_4;
    gamma_3 = params.gamma_3;
    gamma_e = params.gamma_e;

    a = A.^2;
    
    VE = (alpha_1*a)./(beta_1^2 + a).*(1-E) - gamma_1*E.*I - gamma_e*E;
    VA = alpha_2 + alpha_3*E - gamma_2*A;
    VI = gamma_3*(-I + alpha_4*E);  

    
end