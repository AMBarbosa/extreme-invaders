# Extreme invaders: data and scripts

This is the companion repository for the article:<br><br>
<big>Baquero R.A., Barbosa A.M., Ayllón D., Guerra C., Sánchez E., Araújo M.B. & Nicola G.G. (2021) <b>Potential distributions of invasive vertebrates in the Iberian Peninsula under projected changes in climate extreme events</b>. <i>Diversity and Distributions</i>, 27(11): 2262-2276 (https://doi.org/10.1111/ddi.13401)</big>
<br><br>

The environmental variables were also published in a companion data paper:<br><br>
Baquero R.A., Barbosa A.M., Ayllón D., Guerra C., Sánchez E., Araújo M.B. & Nicola G.G. (2021) Environmental data from: Potential distributions of invasive vertebrates in the Iberian Peninsula under projected changes in climate extreme events. Dryad, DOI https://doi.org/10.5061/dryad.rxwdbrv8x
<br><br>

We modelled the distributions of six species based on five regional climate models and their projections for a set of variables representing extreme conditions. We used an ensemble of four modelling methods: Generalized Linear Models (GLM), Generalized Additive Models (GAM), Random Forests, and Bayesian Additive Regression Trees (BART). We assessed model performance using spatial block cross-validation and metrics focused on discrimination and calibration. We finally extrapolated the well-performing models to future climate scenarios, and then quantified and mapped the changes.

The repository contains:
- the values of the regional climatic extreme <b>variables</b> computed by members of the team;
- the redistributable parts of the <b>species occurrence data</b>, namely those currently available on GBIF or published in data papers (citations included in the script and in the article);
- the <b>R scripts</b> which reproduce our complete analysis and resulting plots. The results won't be exactly the same if you run them, as in the article we used several additional data on species occurrence which need to be obtained directly from their original sources upon acceptance of their conditions.
