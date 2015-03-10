% Usage for svmpredict and svmtrain in libsvm
% =====
% 
% matlab> model = svmtrain(training_label_vector, training_instance_matrix [, 'libsvm_options']);
% 
%         -training_label_vector:
%             An m by 1 vector of training labels (type must be double).
%         -training_instance_matrix:
%             An m by n matrix of m training instances with n features.
%             It can be dense or sparse (type must be double).
%         -libsvm_options:
%             A string of training options in the same format as that of LIBSVM.
% 
% matlab> [predicted_label, accuracy, decision_values/prob_estimates] = svmpredict(testing_label_vector, testing_instance_matrix, model [, 'libsvm_options']);
% matlab> [predicted_label] = svmpredict(testing_label_vector, testing_instance_matrix, model [, 'libsvm_options']);
% 
%         -testing_label_vector:
%             An m by 1 vector of prediction labels. If labels of test
%             data are unknown, simply use any random values. (type must be double)
%         -testing_instance_matrix:
%             An m by n matrix of m testing instances with n features.
%             It can be dense or sparse. (type must be double)
%         -model:
%             The output of svmtrain.
%         -libsvm_options:
%             A string of testing options in the same format as that of LIBSVM.
% 
% Returned Model Structure
% ========================
% 
% The 'svmtrain' function returns a model which can be used for future
% prediction.  It is a structure and is organized as [Parameters, nr_class,
% totalSV, rho, Label, ProbA, ProbB, nSV, sv_coef, SVs]:
% 
%         -Parameters: parameters
%         -nr_class: number of classes; = 2 for regression/one-class svm
%         -totalSV: total #SV
%         -rho: -b of the decision function(s) wx+b
%         -Label: label of each class; empty for regression/one-class SVM
%         -ProbA: pairwise probability information; empty if -b 0 or in one-class SVM
%         -ProbB: pairwise probability information; empty if -b 0 or in one-class SVM
%         -nSV: number of SVs for each class; empty for regression/one-class SVM
%         -sv_coef: coefficients for SVs in decision functions
%         -SVs: support vectors
% 
% If you do not use the option '-b 1', ProbA and ProbB are empty
% matrices. If the '-v' option is specified, cross validation is
% conducted and the returned model is just a scalar: cross-validation
% accuracy for classification and mean-squared error for regression.
% 
% More details about this model can be found in LIBSVM FAQ
% (http://www.csie.ntu.edu.tw/~cjlin/libsvm/faq.html) and LIBSVM
% implementation document
% (http://www.csie.ntu.edu.tw/~cjlin/papers/libsvm.pdf).
% 
% Result of Prediction
% ====================
% 
% The function 'svmpredict' has three outputs. The first one,
% predictd_label, is a vector of predicted labels. The second output,
% accuracy, is a vector including accuracy (for classification), mean
% squared error, and squared correlation coefficient (for regression).
% The third is a matrix containing decision values or probability
% estimates (if '-b 1' is specified). If k is the number of classes
% in training data, for decision values, each row includes results of 
% predicting k(k-1)/2 binary-class SVMs. For classification, k = 1 is a
% special case. Decision value +1 is returned for each testing instance,
% instead of an empty vector. For probabilities, each row contains k values
% indicating the probability that the testing instance is in each class.
% Note that the order of classes here is the same as 'Label' field
% in the model structure.