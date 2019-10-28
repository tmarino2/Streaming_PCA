# MB-MSG
Code base is taken from "Stochastic PCA with ℓ2 and ℓ1 Regularization".
Create 4 folders: code, data, page, plot
Download the source code into the code folder, and your data into the data folder.
Run the Demo.m file, which creates a synthetic dataset and runs various stochastic PCA algorithms.

To reproduce the experiments for the synthetic dataset run vrmsgWrapper(k,numiters,dataname,d,gap,num_points) with k=1,3,7 (or any value of choice), numiters=1, dataname = 'synthetic', gap=0.1, num_points = 4*1e4.

To reproduce the experiments for MNIST dataset run vrmsgWrapper with k=1,3,7, numiters=10, dataname='MNIST'.

To run on other datasets, follow the instructions for MNIST, but replace with the dataset name and have the dataset in the "/data" folder.

If you would like to generate synthetic data with other k,gap or d, uncomment syn_gen_exp_gap function in vrmsgWrapper. For more details see the description of the function.
