# Function Descritions

Here are short descriptions of this package functions

In main.R: 
- `GetTRTpred`: Get TRT prediction from best model
- `GetClonePred`: Get the clone-wise score
- `GetModelPrediction`:  Get the TRT prediction from any model
- `TrainModel`: Train any model

In signature_method.R
- `GetPredictionSignature`: Get the signature model prediction when given the model
- `GetSignature`: Get signature from DE results
- `GetSignatureList`: Get signature.list from DE results
- `SignatureCrossValidation`: Perform the inner-workings of the Cross-Validation for the signature scoring model

In Logistic_regression.R
- `GetLRPrediction`: Get Logistic Regression predictions
- `LRCrossValidation`: Perform the inner-workings of the Cross-Validation for the LR model

In Data_preparation.R
- `PrepareData`: Overall function to prepare data for TRTpred pipeline
- `FeatureTrans`: Transform the feature into another feature space
- `DiscrAnalysisBinary`: Binary discriminant analysis function
- `RemoveMutuallyCorrFeatures`: Remove mutually correlated features in matrix
- `PrepareTrainingDataFromSeurat`: Prepare data from Seurat object

In Cross_Validation.R
- `CreateCVFolds`: Create Cross-Validation (CV) folds
- `CreateNestedCVFolds`:  Create Nested-Cross-Validation (NCV) folds
- `GetBestHyperparams`: Get the best hyperparameters from any cross-validataion 
- `NestedCrossValidation`: Run the Nested-Cross-Validation (NCV)
- `CrossValidation`:  Run the Nested-Cross-Validation (CV)

In DE_function.R
- `getPseudoBulk`: Function to get the pseudo-bulk matrix
- `RunDEA`: Function to run Differential Expression analysis
- `Runlimma`: Run limma (voom/trend) Differential Expression analysis
- `RunedgeR`: Run edgeR (LRT/QFL) Differential Expression analysis
- `RunDESeq2`:Run DESeq2 (LRT/QFL) Differential Expression analysis
- `RunSeurat`:  Run Seurat (LR/Wilcox/NegBin) Differential Expression analysis (findMarkers)

In Model_evaluation
- `BinarizeScorePrediction`: Binarize the continuous score into TRT predictions (TRUE / FALSE)
- `GetAccuracyMetrics`: Get the overal metrics (the one in the functions below)
- `GetMetricsBinary`: Get the TP, TN, FP, FN
- `GetMccBinary`: Get the MCC
- `GetAccurayBinary`: Get the Accuracy
- `GetF1Binary`: Get the F1
- `GetKappaBinary`: Get the Kappa
- `GetAUC`: Get the AUC
- `GetSensitivityBinary`: Get the sensitivity 
- `GetSpecificityBinary`: Get the specificity

In Signature_score_function.R
- `GetSignatureScore`: Get signature score
- `RunSigScoreAverage`: Run signature score average
- `RunSigScoreUCell`: Run signature score UCell
- `RunSigScoreAUCell`: Run signature score AUCell
- `RunSigScoreSingscore`: Run signature score Singscore
- `RunSigScoreScGSEA`: Run signature score scGSEA
- `run_scgsea`: run inner workings of scGSEA

In helpers.R
- `ProcessConfigFile`: Process configuration files describing the models
