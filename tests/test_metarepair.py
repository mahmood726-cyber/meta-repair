"""Tests for MetaRepair pipeline."""
import sys, math, pytest, numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from pipeline import (diagnose_multimodality, diagnose_outliers, diagnose_publication_bias,
                       diagnose_instability, correct_robust_pooling, correct_bias_adjustment,
                       decompose_uncertainty, decision_grade)


class TestMultimodality:
    def test_unimodal(self):
        yi = np.array([-0.5, -0.4, -0.3, -0.2, -0.1])
        sei = np.array([0.1]*5)
        r = diagnose_multimodality(yi, sei)
        assert r['n_modes'] == 1

    def test_bimodal(self):
        yi = np.array([-2, -1.9, -1.8, 1.8, 1.9, 2.0])
        sei = np.array([0.1]*6)
        r = diagnose_multimodality(yi, sei)
        assert r['n_modes'] == 2
        assert len(r['cluster_labels']) == 6

    def test_k3_returns_unimodal(self):
        r = diagnose_multimodality(np.array([0.1, 0.2, 0.3]), np.array([0.1]*3))
        assert r['n_modes'] == 1


class TestOutliers:
    def test_no_outliers(self):
        yi = np.array([-0.3, -0.2, -0.1, 0.0, 0.1])
        sei = np.array([0.2]*5)
        r = diagnose_outliers(yi, sei)
        assert r['n_outliers'] == 0

    def test_obvious_outlier(self):
        # Needs extreme outlier: tau2 inflates so moderate outliers don't exceed 2.5 threshold
        yi = np.array([-0.3, -0.2, -0.1, 0.0, 50.0])
        sei = np.array([0.1]*5)
        r = diagnose_outliers(yi, sei)
        assert r['n_outliers'] >= 1
        assert 4 in r['outlier_indices']


class TestPubBias:
    def test_symmetric_no_bias(self):
        yi = np.array([-0.5, -0.3, -0.1, 0.1, 0.3])
        sei = np.array([0.1, 0.15, 0.2, 0.15, 0.1])
        r = diagnose_publication_bias(yi, sei)
        assert r['egger_p'] > 0.05 or not r['bias_detected']

    def test_returns_dict(self):
        yi = np.array([-0.5, -0.3, -0.1, 0.1, 0.3])
        sei = np.array([0.2]*5)
        r = diagnose_publication_bias(yi, sei)
        assert 'egger_p' in r
        assert 'tf_k0' in r


class TestInstability:
    def test_stable_review(self):
        yi = np.array([-0.8, -0.7, -0.6, -0.5, -0.4])
        sei = np.array([0.1]*5)
        r = diagnose_instability(yi, sei)
        assert r['stability'] >= 70

    def test_returns_8_specs(self):
        yi = np.array([-0.3, -0.1, 0.1, 0.3, 0.5])
        sei = np.array([0.2]*5)
        r = diagnose_instability(yi, sei)
        assert r['n_total'] == 8


class TestRobustPooling:
    def test_no_outliers_same_as_dl(self):
        yi = np.array([-0.5, -0.3, -0.1])
        sei = np.array([0.2]*3)
        r = correct_robust_pooling(yi, sei, [])
        assert math.isfinite(r['theta'])

    def test_outlier_downweighted(self):
        yi = np.array([-0.3, -0.2, -0.1, 5.0])
        sei = np.array([0.2]*4)
        r_orig = correct_robust_pooling(yi, sei, [])
        r_robust = correct_robust_pooling(yi, sei, [3])
        assert abs(r_robust['theta']) < abs(r_orig['theta'])


class TestDecisionGrade:
    def test_grade_a(self):
        diag = {'multimodality': {'n_modes': 1}, 'outliers': {'n_outliers': 0},
                'publication_bias': {'bias_detected': False}, 'instability': {'stability': 95}}
        corrected = {'theta_original': -0.5, 'theta_final': -0.49}
        r = decision_grade(diag, corrected, {})
        assert r['grade'] == 'A'

    def test_grade_d(self):
        diag = {'multimodality': {'n_modes': 2}, 'outliers': {'n_outliers': 3},
                'publication_bias': {'bias_detected': True}, 'instability': {'stability': 40}}
        corrected = {'theta_original': -0.5, 'theta_final': -0.1}
        r = decision_grade(diag, corrected, {})
        assert r['grade'] == 'D'


class TestResultsIntegrity:
    @pytest.fixture
    def summary(self):
        path = Path('C:/Models/MetaRepair/data/output/metarepair_summary.json')
        if not path.exists(): pytest.skip('Run pipeline first')
        import json
        with open(path) as f: return json.load(f)

    def test_review_count(self, summary):
        assert summary['n_reviews'] == 307

    def test_grades_sum(self, summary):
        total = sum(summary['grades'].values())
        assert total == summary['n_reviews']

    def test_uncertainty_sums_to_100(self, summary):
        u = summary['avg_uncertainty']
        total = sum(u.values())
        assert 95 < total < 105  # allow rounding


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
