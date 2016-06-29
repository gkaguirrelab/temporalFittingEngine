function [A,Ainv,Q] = PoirsonWandellEllipsoid(conditionStr,sf)
% [A,Ainv,Q] = PoirsonWandellEllipsoid(conidtionStr)
%
% Return the parameters in cone contrast space of the threshold ellipsoid
% predicted by the Poirson-Wandell model.
%    Poirson AB, Wandell BA. 1996. Pattern-color separable pathways predict
%    sensitivity to simple colored patterns Vision Res 36: 515-26.
%
% 6/29/16  dhb  Wrote it.

%% Get parameters from paper Tables 1 and 2
%
% These give the model parameters for each of the three subject/condition
% combinations measured.
switch (conditionStr)
    case 'HT,cc'
        % Subject HT, constant cycles stimulus as sf varies
        MConesToOpponent = [0.759 0.649 0.058 ; -0.653 0.756 0.033 ; -0.159 -0.414 0.896];
        theSfs = [0.5 1 2 4 ]';
        theFactors = [0.753 8.732 4.075 ;
            1.226 8.044 2.567 ;
            0.864 5.435 1.187 ;
            0.683 2.536 0.266];
    case 'HT,cs'
        % Subject HT, constant size stimulus as sf varies
        MConesToOpponent = [-0.246 0.967 -0.072 ; -0.698 0.716 -0.005 ; 0.165 -0.599 -0.784];
        theSfs = [0.5 1 2 4 8]';
        theFactors = [11.864 63.945 8.080 ;
            18.037 74.214 7.212 ;
            23.739 60.598 3.235 ;
            28.354 31.746 3.379 ;
            26.959 12.302 2.622];
    otherwise
        error('Unknown condition string entered')
end

index = find(theSfs == sf);
if (isempty(index))
    error('No data at specified spatial frequency');
end
    
A = diag(theFactors(index,:))*MConesToOpponent;
Ainv = inv(A);
Q = A'*A;