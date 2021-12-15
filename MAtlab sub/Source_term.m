function[f] = Source_term(f_term,J)
%% create the local element source term vector to be called into the global vector 
f = f_term*J*[1;1];
%can also be achieved with ones[2:1]*f_term*J
end