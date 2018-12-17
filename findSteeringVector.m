function [IndexAt,IndexAr] = findSteeringVector(H,At,Ar,Ns)
%%% This is the function to obtain the beam steering vectors at the transmitter and receiver 
%%% by finding the array response vectors corresponding to the largest effective channel gain (an exhaustive search).
%%% Only suitable for single data stream case.
% 
% Input
%       H     : channel matrix, Nr x Nt
%       At    : the collection of transmit steering vectors, Nt x NumPath
%       Ar    : the collection of receive steering vectors, Nr x NumPath
% Output
%      IndexAt : index of the path selected at the transmitter
%      IndexAr : index of the path selected at the receiver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nt = size(At,1);
Nr = size(Ar,1);

NumPath = size(At,2);

EffGain = zeros(NumPath);

if Ns == 1    
    
    for t = 1:NumPath
        
        for r = 1:NumPath
        
            EffGain(t,r) = abs(Ar(:,r)'*H*At(:,t));
        
        end
    
    end
else
    
    error('Not single data stream.');
    
end


[row,col] = find(EffGain == max(max(EffGain)));

IndexAt = row;
IndexAr = col;
      
end


