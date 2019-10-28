%% Generate and save a synthetic data set
% k is the component at which the gap is
% d is the dimensionality of the data
% gap is the size of the gap
% num_points is the size of the data set
function syn_gen_exp_gap( k,d,gap, num_points )
    mu = zeros(1,d);
    step = .1;
    Sigma = ones(1,d);
    Sigma(k+1:end) = gap*2.^(0:-step:-(d-k-1)*step);
    X = mvnrnd(mu,Sigma,num_points);
    X = X';
    data.training=X(:,1:.8*num_points);
    data.tuning=X(:,.8*num_points+1:.9*num_points);
    data.testing=X(:,.9*num_points+1:num_points);
    fname = sprintf('../data/syn_k=%d_d=%d_gap=%0.2e_num_points=%d.mat',k,d,gap,num_points);
    save(fname,'data');
end

