function [ x ] = create_piecewise_signal( t_sig, itk, a_k )
%CREATE_PIECEWISE_SIGNAL Create a piecewise constant signal based on the
%discontinuities located at itk, of amplitudes a_k

no_pieces = length(itk);
for i = 2:length(a_k)
    a_k(i) = a_k(i)+ a_k(i-1);
end
% a_k(2:end) = a_k(1:end-1)+a_k(2:end);

% Generate the continuous-time signal x(t)
x = zeros(size(t_sig));
for i = 1:no_pieces-1
    x(itk(1,i): itk(1,i+1)) = a_k(1,i);
end
x(itk(1,no_pieces):end) = a_k(1,no_pieces);

if no_pieces == 1
    x(itk: end) = a_k;
end



end

