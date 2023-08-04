# LassaFever_SeasonalityNigeria

Welcome to the GitHub repository for the research paper _Modelling seasonality of Lassa fever in Nigeria_. This repository contains the MATLAB code and data used to produce the results and figures presented in the paper. Below is an overview of the repository structure and instructions on how to navigate and use the contents.

## Repository Structure

The repository is organized into several folders, each serving a specific purpose:

1. **Fits_Models**: Contains the MATLAB functions used to run the model developed in the paper and perform model fitting using an Approximate Bayesian Computation (ABC) scheme. It includes the core codebase necessary for reproducing the results in the paper.

2. **GraphScripts**: In this folder, you'll find the MATLAB scripts that were employed to generate the figures showcased in the research paper. These scripts utilise the results obtained from the model fitting process.

3. **Inputs**: The **Inputs** folder holds the input data and configurations required to fit the model using ABC algorithms. Data_TimeSeries.mat is the data as a time series that was used to fit the data. 5 and 6 parameter inputs include parameter priors, initial conditions, and other input files necessary to run the fitting process.

4. **Results**: The **Results** folder contains the outputs generated from the ABC fitting procedure. This includes the posterior samples, summary statistics, and any other relevant data produced during the model fitting. These results are then used in the graph generation process.

## Citation

If you find this work helpful or build upon it, please consider citing our original paper:

\[J. McKendrick, W. Tennant, M. J. Tildesley. "Modeling seasonality of Lassa fever in Nigeria"\]

## Contact

For any inquiries or assistance, you can contact us at j.mckendrick@warwick.ac.uk.
