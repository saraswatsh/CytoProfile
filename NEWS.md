# CytoProfile 0.1.2

* Added references in description field of the DESCRIPTION file to the multivariate methods mentioned in the format 'Authors (Year) <doi:...>'.

* 'cyt_rf', 'cyt_splsda', 'cyt_ttest', 'cyt_dualflashplot', and 'cyt_volc' now have a logical 'verbose' argument to print the output if user wants to. This is to ensure printed output in console can be easily suppressed.

* 'cyt_xgb' has a new logical argument called 'print_results' that works the same as the 'verbose' argument in other functions mentioned above. 

* Have removed the arguments for changing graphical parameters in 'cyt_bp', and 'cyt_bp2' so it no longer changes the user's graphical parameters. 

* Updated examples from graphical parameters changing without resetting properly to now the graphical parameters reverting to original using 'on.exit()' after the execution of code.

* Updated 'getting_started.Rmd' vignette where the files created and saved are now saved to a temporary directory using 'tempdir()' to avoid creating files in the user's working directory. Additionally, 
the vignette now uses 'on.exit()' to revert to original graphical parameters for the 'cyt_errbp' examples shown.

* 'cyt_rf', 'cyt_splsda', and 'cyt_xgb' now has a 'seed' argument to set the seed for reproducibility. 

# CytoProfile 0.1.1

* Fixed some general CRAN automated test issues. 

* Updated functions and examples to avoid generating PDF and PNG files during 
checks.

* Functions now also generated plots within RStudio graphics device and provides 
the option to save the plots as PDF or PNG files depending on the function.

# CytoProfile 0.1.0 

* Initial CRAN submission.

Since this is our initial submission, our NEWS.md contains a single entry 
indicating the initial release. We plan to provide more detailed change 
logs in subsequent versions.

# CytoProfile 0.0.0.9000 (Development Version) 

* Preparation For Initial CRAN Submission

This version is the development version of CytoProfile. We are preparing
for the initial CRAN submission. Below are the features added in this version:

## Major Changes

- Added functions for ANOVA, Boxplots, Dual-Flashlight Plot, Error-Bar Plot,
Heatmap, Principal Component Analysis (PCA), Random Forest, Skewness and Kurtosis, 
Sparse Partial Least Squares - Discriminant Analysis (sPLS-DA), T-Test, 
Volcano Plot, and Extreme Gradient Boosting (XGBoost).

## Minor Changes and Bug Fixes

- Added a vignette to demonstrate the usage of the package.
- Added built-in data sets for demonstration purposes.
- Examples from the documentation are now working as expected.
- 'cyt_splsda' is now able to use one grouping column to run analysis. 
- 'cyt_skku' now prints the output in readable format.
