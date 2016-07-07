function ellParams = QToEllParams(Q)
% ellParams = QToEllParams(Q)
%
% Take a positive-definite symmatrix matrix Q and convert it to the six
% parameters of an ellipsoid, in the ellParams form of three scalars and
% three euler angles, as one column vector.
%
% Notice that the may recover the scalar prameters in a different order than
% we put them in with, so that creating Q from ellParams and then recovering
% the ellParams vector can lead to two different parameter vectors.  What is 
% preserved is that both paramter vectors lead to the same Q. 
%
% This ambiguity doesn't both us, in part because it is not clear we really
% need this routine for any serious purpose -- the parameter vectors are
% just for the search routines, and we don't really ever need to get them
% back once we have Q.
%
% 7/4/16  dhb  Wrote it.

[u,s,v] = svd(Q);
scalers = 1./diag(sqrt(s))';
eul = rotm2eul(u);
ellParams = [scalers  eul]';
[ACheck,AinvCheck,QCheck] = EllipsoidMatricesGenerate(ellParams);
if (max(abs(QCheck(:)-Q(:))) > 1e-8)
    error('Cannot recover ellipsoid parameters from Q');
end