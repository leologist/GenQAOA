function psi = getProductState(angles)
%getProductState returns a product state of qubits, individually oriented
%                on the bloch sphere with (theta, phi) angles
%                specified by the rows of the input

    N = size(angles,1);
    bloch = @(theta, phi) [cos(theta/2); exp(1i*phi)*sin(theta/2)];

    theta = angles(1, 1);    phi = angles(1, 2);
    psi = bloch(theta, phi); 
    for i = 2:N
        theta = angles(i, 1);   phi = angles(i, 2);
        psi = kron(psi, bloch(theta,phi));
    end
 
end