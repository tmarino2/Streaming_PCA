
%% Extracts the best maximal kk-dimensional subspace from a solution
% found by pca_update_F or pca_update_vN
%
% kk - the dimension of the subspace which we seek
% vectors, values - "nontrivial" eigenvectors and eigenvalues
%                   of the iterate
%%
function best_vectors = pca_solution_original( kk, vectors, values )
if ( length( values ) ~= size( vectors, 2 ) )
    error(['pca_solution : number of columns of vectors must',...
        ' be same as size of values']);
end
if ( length( values ) < kk )
    error(['pca_solution : must have at least a kk-dimensional',...
        ' nontrivial subspace']);
end
[ ~, indices ] = sort( values, 'descend' );
best_vectors = vectors( :, indices( 1:kk ) );
end
