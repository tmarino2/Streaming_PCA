for i = [1 3 7]
            k = i;
            numiters = 10;
            fname = sprintf('../data/MNIST.mat');
            load(fname);
            
            methods={'msg','l2rmsg','sgd','vrmsg','l2vrmsg'};
            method_names = {'MSG','RMSG','Oja','VR-MSG','VR-RMSG'};
            eta = .5;
            batch_const = 10;
            RST = 1;

            N=size(data.training,2);
            d = size(data.training,1);
            num_points = N;
            data.training=data.training(:,1:N);
            mu=mean(data.training,2);
            data.training=data.training-repmat(mu,1,N);
            data.tuning=data.tuning-repmat(mu,1,size(data.tuning,2));
            data.testing=data.testing-repmat(mu,1,size(data.testing,2));
            lambda_k = svds(data.training/sqrt(N),k);
            data.training = data.training/lambda_k(k);


            S = svds(data.training/sqrt(N),k+1).^2;
            S = sort(S,'descend');
            gap = S(k) - S(k+1); 
            beta = (S(k) + S(k+1))/4;
            lambda = gap/2;
            
            dataname = 'MNIST';
            plotobjV(dataname,methods,k,N,numiters,eta,lambda,beta,method_names);
end