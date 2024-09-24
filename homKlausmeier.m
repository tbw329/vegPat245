function x=homKlausmeier(W,B,P,L,R,J,M) 
x(1) = P - L*W - R*W*B.^2;
x(2) = J*R*W*B.^2 - M*B.^2;
end
