# Addressing prior feedback
- We took Fred's feedback and appended our proposal to clarify software we will be using (ggplot, dplyr, geomorph) and clarified how the work will be split up. Rashad will organize the data and make PCA plots while Sara will run statistical analysis on our GPA across different groups. There will be a lot to dig into there.

# Progress since last submission
We have produced an R script that does the following tasks:
- Adds beetle data from disparate files into two master TPS files - Elytron.TPS and Pronotum.TPS.
- When added to the master file, beetles retain the name of the original file they came from, allowing us to identify their tribe.
- Each individual beetle also has a new number when added to the master file so each one can be identified with a unique number.
- Using the master files, we performed generalized procrustes analysis on each to eliminate size as a variable. This allows us to analyze shape alone.
- We plotted the pronota and elytra on two PCA plots shown below.
![Hopefully the figure shows up here :)](../Project_work/Pronotum/Pronotum_plot.png)
![Hopefully the figure shows up here too :)](../Project_work/Elytron/Elytron_plot.png)

# Future steps
- With GPAs run and PCAs made we can run statistical tests like ANOVA to detect significant differences in shape across different groups of beetles. We can also cross-reference our results with a phylogeny to see what patterns or interesting data emerge.
- We also want to be able to plot the landmarks from pronota and elytra of individual beetles to get an idea of what changes each PC is visually describing. Looking into interactive plot packages in R for this.

# Project organization
Repository contains a handful of files:
- Data - contains the original data from the paper including zipped and unzipped files of landmark and classification info.
- OG_Paper - contains the original paper as a PDF as well as figures from the paper that we reference in our README.md.
- Project_work - contains our code, master TPS files, and figures, separated into files for the elytra and pronota.
- doc - contains our check-ins.
- README.md 

# Struggles
- Initially struggled to figure out how to run GPA and make PCAs with all the data the way the authors had it. We fixed this by consolidating them into two master files where each sample retains the identity of its original file.
- We struggled with how we could pull out the unique numerical identifiers for each beetle once they were added. This was solved by adding a # before each one that we could then grep for.
- Bugs were encountered and fixed here and there throughout the whole process of coding.