import numpy as np
from sklearn.base import TransformerMixin
from sklearn import mixture


class BiBiTransformer(TransformerMixin):
    """Obtain binormality/bimodality information for a given a distribution.
    This is accomplished by fitting the dataset to two Gaussian Mixture models with one and two components, respectively.
    Then the BIC score together with Bayes Factor is used to assert binormality, and a modification of Helguero's theorem to assert bimodality.
    Other parameters such as means, variances and weights are given for both distributions.

    :param float confidence_level: confidence level to use in Bayes Factor when determining binormality. (Default 5.0)
    :param **kwargs: other arguments to be passed to sklearn.mixture.GaussianMixture

    :raises ValueError: if ``confidence_level`` is not between 0 and 100.
    """

    def __init__(self, confidence_level=5.0, **kwargs):
        self.models = self.GMM_models(**kwargs)

        if confidence_level < 0 and confidence_level > 100:
            raise ValueError("confidence_level must be between 0 and 100")
        self.confidence_level = confidence_level

    def GMM_models(self, **kwargs):
        """Create two GaussianMixture models with the same hyperparameters,
        one with one component and the other with two components.

        :return: tuple of both Gaussian Mixture models
        """
        gmm_1 = mixture.GaussianMixture(
            n_components=1,
            **kwargs)
        gmm_2 = mixture.GaussianMixture(
            n_components=2,
            **kwargs)
        return gmm_1, gmm_2

    def fit(self, X):
        """Fits both models to the same data.

        :param X: array of shape (n_samples, n_features)"""
        self.X = X
        self.models = tuple(map(
            lambda model: model.fit(X),
            self.models))
        return self

    def transform(self, X, y=None):
        """Get parameters describing the distribution fitted to model with lowest BIC value.
            The information returned is: 
            bimodality
            binormality/uninormality/insuficient evidence
            BIC values
            mean(s)
            variance(s)
            weight(s)

        :param X: array of shape (n_samples, n_features)

        :return: Dictionary with distribution information."""
        bics = []
        means = []
        variances = []
        weights = []
        for gmm in self.models:
            m = gmm.means_.flatten()
            v = gmm.covariances_.flatten()
            b = gmm.bic(X)
            w = gmm.weights_.flatten()
            means.append(m)
            variances.append(v)
            bics.append(b)
            weights.append(w)

        # bayes factor criteria for normality
        uninormal, binormal, insuf_ev, p = self.bayes_factor_criteria(bics[0],bics[1])
        if binormal:
            maxm = np.argmax(means[1])
            minm = np.argmin(means[1])
            mean1 = means[1][minm]
            var1 = variances[1][minm]
            w1 = weights[1][minm]
            mean2 = means[1][maxm]
            var2 = variances[1][maxm]
            w2 = weights[1][maxm]
            # Helgueros theorem for unimodality
            unimodal = self.is_unimodal(
                [mean1, mean2], [var1, var2])
        else:
            mean1 = means[0][0]
            var1 = variances[0][0]
            w1 = weights[0][0]
            mean2, var2, w2 = np.nan, np.nan, 0
            unimodal = True
        info = dict(
            binormal=binormal,
            uninormal=uninormal,
            insuf_ev=insuf_ev,
            unimodal=unimodal,
            bics=bics,
            p=p,
            mean1=mean1,
            mean2=mean2,
            var1=var1,
            var2=var2,
            w1=w1,
            w2=w2)
        return info

    def bayes_factor_criteria(self, bics0, bics1):
        """
        Determine uninormality, bimodality and insufficient evidence parameters from Bayes Factor.

        :param bics0: BIC value for 1-component model
        :param bics1: BIC value for 2-component model

        :return: tuple with uninormal, binormal and insuficient evidence boolean values.
        """
        diff_bic = bics1 - bics0
        # probability of a two-component model
        p = 1 / (1+np.exp(0.5*diff_bic))
        if p == np.nan:
            if bics0 == np.nan:
                p = 1
            elif bics1 == np.nan:
                p = 0

        uninormal = p < (self.confidence_level / 100)
        binormal = p > (1 - (self.confidence_level / 100))
        insuf_ev = True if (not uninormal and not binormal) else False
        return uninormal, binormal, insuf_ev, p

    def is_unimodal(self, means, vars):
        """implements P.Dans version of Helguero's theorem in order to detect unimodality.

        :param means: array with values of means for both fitted models.
        :param vars: array with values of variances for both fitted models.

        :return: unimodal boolean value."""
        r = vars[0] / vars[1]

        # separation factor
        s = np.sqrt(
            -2 + 3*r + 3*r**2 - 2*r**3 + 2*(1 - r + r**2)**1.5
        ) / (
            np.sqrt(r)*(1+np.sqrt(r))
        )

        unimodal = abs(means[1]-means[0]) <= s * \
            (np.sqrt(vars[0]) + np.sqrt(vars[1]))
        return unimodal
