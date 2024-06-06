# seqwrap - Item-by-item Iterative Model Fitting


The `seqwrap` R package allows you to model high dimensional data using
an item-by-item strategy with a specific model engine. Using `seqwrap`
you can fit e.g. gene/transcript data from a RNA sequencing experiment
to a model specified with random effects or custom distributions. This
is possible as `seqwrap` efficiently iterates over all item (e.g. genes)
and fits them to the same model formulation using any model fitting
algorithm implemented in R. The package also provides pre- and
post-fitting convenience functions. We show how to implement custom
models to commonly analysed data and how pre- and post-fitting
processing can be used to share information across independent models.
