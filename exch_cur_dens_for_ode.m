%% Exchange Current Density function i_0
%   Created July 12, 2011 by Scott Moura, 
%   Modified July 22, 2024 by Yuichi Kajiura 

function [i_0n,i_0p,varargout] = exch_cur_dens_for_ode(p,c_ss_n,c_ss_p,c_e)

% Compute exchange current density
i_0n = p.k_n * ((p.c_s_n_max - c_ss_n) .* c_ss_n .* c_e).^p.alph;
i_0p = p.k_p * ((p.c_s_p_max - c_ss_p) .* c_ss_p .* c_e).^p.alph;
di_0n = p.k_n * (p.alph*((p.c_s_n_max - c_ss_n) .* c_ss_n .* c_e).^(p.alph-1).*( c_e.*(p.c_s_n_max - c_ss_n) -  c_ss_n .* c_e) );
di_0p = p.k_p * (p.alph*(max((p.c_s_p_max - c_ss_p),0) .* c_ss_p .* c_e).^(p.alph-1).*( c_e.*(p.c_s_p_max - c_ss_p) -  c_ss_p .* c_e) );

if(nargout > 2)
    varargout{1}=di_0n;
    varargout{2}=di_0p;
end