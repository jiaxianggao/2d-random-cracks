function internalV = calculateInternalV(X)
    Node1 = evalin('base','Node1');
    Element1=evalin('base','Element1');
    I=evalin('base','I');
    z=evalin('base','z');
    MeasPatt=evalin('base','MeasPatt');
    cols = size(X,2);
    internalV = zeros(5383,108,cols);
    for col = 1:cols
        sigma_data=X(:,col);
        vh = ForwardSolution(Node1,Element1,I,sigma_data,z,MeasPatt,'real');
        internalV(:,:,col)=vh.Current; 
    end    
end