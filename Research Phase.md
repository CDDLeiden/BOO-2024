# BOO 2024 Research Phase

The idea of this research phase is to simulate the basics of bioactivity modeling in a short time span. Using the in-house 
**[QSPRpred](https://github.com/CDDLeiden/QSPRpred)** package, we can create a regression model that can predict the bioactivity
for a particular protein and subsequently subject it to a virtual screen to visualise if there are any molecules that 
perform well on the protein that we have not seen tested in literature before. An extra bonus step would be to compare 
the different proteins in a family to see how much specificity is found in the compounds. An example of this process can 
be found in **[this article](https://pubs.acs.org/doi/10.1021/acs.jcim.2c01645)** for example.

This page is intended as a small guide to the big steps in this project, to give you an idea of what is asked of you, 
and the background behind it. It will be more of a guide to keep you on track, than an actual step-by-step tutorial.

## Step 0: Python
To interface with QSPRpred and related packages, knowledge of Python is required. There are plenty of resources available 
to learn Python, as it is arguably one of the easiest coding languages to pick up. Coding does need a different way of 
looking at problems, and it might take some time to get into that mindset - don't be discouraged! You can learn a lot by 
just playing around with your code, and often the runtime of your code can be short enough to see what effects certain 
changes had made. It is encouraged to not try to follow all the code examples to the letter, but to understand what functions 
or variables mean.

Aside from the provided **[Python Tutorial](https://github.com/jesperswillem/python-course-materials)** from Leiden University, 
the following resources can help out as well for some extra training / reading.

- [Jake VanderPlas's notebooks on Data Science](https://github.com/jakevdp/PythonDataScienceHandbook) - Focuses specifically
 on some of the most important packages such as NumPy, Pandas and Scikit-Learn
- [IPython Cookbook](https://ipython-books.github.io/) - Shows how to perform data science on Jupyter Notebook
- [Chris Albon's notes on machine learning / AI](https://chrisalbon.com/Home) - Collected notes and insights to help you
 a different perspective
- [Vincent D. Warmerdam's Scikit-learn Crash Course](https://www.youtube.com/watch?v=0B5eIE_1vpU) - If you are more of a 
 visual learner, Vincent has made an excellent, comprehensive tutorial on how to perform machine learning using Scikit-Learn.

## Step 1: Creating a Regression Model
The overarching goal is to predict the bioactivity for new compounds on a specified protein. However, to do that, a model 
that can predict these bio-activities is needed. In essence, this comes down to the following three major steps:
- Curate a training dataset that contains an **X** (labels, descriptors) and an **Y** (bioactivity value) .
- Create a machine learning model that is trained on the previous dataset.
- Using a test set, validate the machine learning model to show its prediction strength.

Here is a brief overview of each of the steps.

### Training dataset
An often underestimated step of the ML pipeline is the curation and preparing of data. In this project, we are exclusively
using bioactivity data, originating from **[Papyrus](https://github.com/OlivierBeq/Papyrus-scripts)**. The qualitative 
property here is the _pChEMBL_ value, which is an aggregated value for the strength of a compound-protein interaction on a
negative logarithmic scale. An IC50 value of 1 µM (10<sup>-6</sup> mol) will translate to a pChEMBL value of 6. This will be
the aforementioned **Y**. However, this still needs extra data on the properties of the compound (The **X**). This always 
includes a label column (Name of the compound), and several properties that describe the compound in a computer-readable way.
The most common way to do this is to calculate descriptors, such as ECFP6. Combining the **X** and **Y** will provide you
with a compatible dataset for a machine learning model.

To perform a proper (regression) model, a training, validation and test set should be created; often done in an 80/10/10% split;
though this can vary depending on the amount of data available. 
- A training set is needed to let the model learn the relationships that eventually determine what it will predict.
- A validation set is needed to (internally) estimate the prediction error for the model.
- A test set is needed to assess how well the model predicts, as well as check for any overfitting/underfitting. The test set
should only be introduced on the very last step of model creation.

### Machine Learning model
Once the dataset is in place, the actual training of the model can begin. In QSPRpred, setting up a regression model is
relatively simple, though the [tutorials](https://github.com/CDDLeiden/QSPRpred/tree/main/tutorials) present in QSPRpred
go in more detail:
```bash
model = SklearnModel(
    base_dir='<Output directory of your model>',
    alg=KNeighborsRegressor,
    name='RegressionModel',
)
```
There are, however, many places where one can tweak with the parameters (or even optimise them), to gain different results.
On top of that, we chose a K Neighbours Regressor, yet there are plenty of different algorithms available. Experimentation
can lead to better predictive power from your model, and you will find that different algrorithms will work better/worse on
different datasets.

### Validation
Finally, to make sure that the model is up to standards, both internal and external validation can be performed. Internal
validation can be done using Cross-Validation, using the created validation set. The goal here is to assess the prediction
error of the model, to gain a read on how well it is performing internally. For External validation, you can use the test
set. The idea there is to use a set that has not been introduced to the model whatsoever, and predict its bioactivity using
the constructed model. Plotting the predicted values versus the actual values, a plot can be constructed with an associated
R<sup>2</sup> and RMSE value that can show the strength of the predictions.

Once this is all in place, optimisation and tuning of the model (or tweaking with the dataset) can be performed to potentially
increase the predictive power of the model. Once satisfied, it is time to predict molecules that do not have an associated 
bioactivity value yet.

## Step 2: Virtual screening

Model creation is only one half of the equation, and given all the effort needed to create such a model specifically for a
protein, it is only natural that it is put to the test. The two major steps here are:
- Creating a virtual screening dataset and predicting the bioactivity.
- Assessment of high-scoring compounds

### Bioactivity predictions
(Virtual) compound libraries are available online, [Enamine](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7593547/) 
being an example of one that was used in the research mentioned at the beginning. These compound libraries contain a plethora
of compounds; some libraries only pertaining to already synthesised compounds, others to include every possible combination.
Preparing a virtual screening dataset can often take long, as having to calculate descriptors for each compound can stack
if there are many compounds available. Once your dataset is prepared in the same way as your previous dataset, the next step
is to simply predict the bio-activities of the virtual screening compounds.

For this research, a small virtual set has been composed. The virtual screening set is obtained using NextMove Software’s SmallWorld,
a tool for performing similarity-based queries of compounds based on chemical modifications. The obtained dataset consists
of an extensive series of compounds that are analogous to actives for targets investigated in the BOO 2024 projects. The
13470 compounds in the dataset are commercially available through two different chemical vendors, Enamine and eMolecules.

### Quality Assessment
Predicting the bioactivity is only one part of the virtual screen; assessing which compounds to pick for further investigation
is just as important. The simplest way to do this is to rank the predicted compounds based on bioactivity; higher scoring 
compounds likely have more chance to be interesting when tested in an experimental setting. Common practice is to compare
these high scoring compounds to compounds that have already proven to be effective; whether that are high affinity compounds
described in literature, or even drugs that have already hit the market. Assessing the structure and identifying important
functional groups in the compound can lead to more insights in compound-protein interactions. This can be even further
investigated using a 3D technique such as Docking.

An extra approach that can be employed is clustering. If you are ranking compounds based on bioactivity scores only, you
will find that the top is one compound, and many derivatives of that compound. While it makes sense, this is less interesting
when you want to find several different type of compounds, and make your predictions chemically diverse. Clustering the compounds
based on tanimoto distance, using a technique such as T-SNE and picking the centers of these clusters can create a diverse set
of compounds that can unlock even further insights on the compound-protein interaction. 

## Further information

Once that is done, you have created a regression model for a particular protein and screened a virtual set to find novel 
compounds that could potentially work well on that protein. Congratulations! While this is only the base for these types
of studies, it is a good start to understand the concepts behind it - many of the steps that expand on these concepts or 
introduce new ones are related to the ones performed here. 

Listed below are some of the potential follow-up steps:
- 3D techniques: Docking, and even subsequent techniques such as molecular dynamics or free energy pertubation
can give tremendous insights on the interactions found in compound-protein complexes, however these are often dependent 
on the availability of 3D structures.
- Proteochemometrics: This research is performed on only a single protein, however these proteins are often in a family of 
similar(-ish) proteins. Combining information from this entire family (or even more) can often result in an increase of 
both prediction power, and insights gained.
- Deep Learning (not feasible in Noteable): The presented regressor works well and quick, but many newer techniques use the power of reinforcement
training to find deeper relations and increase predictive power; though often at the cost of time and computational power.
