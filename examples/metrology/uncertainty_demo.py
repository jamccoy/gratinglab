"""
Test script to demonstrate Priority 1 improvements
Shows what the new uncertainty tracking provides
"""
import numpy as np

def demonstrate_uncertainty_decomposition():
    """
    Demonstrate how the new uncertainty decomposition works
    """
    
    print("="*70)
    print("PRIORITY 1: FIT UNCERTAINTY TRACKING - DEMONSTRATION")
    print("="*70)
    
    # Simulate some typical AFM blaze angle data
    # Let's say we measured 10 grooves
    true_angles = [17.1, 17.3, 16.9, 17.2, 17.0, 17.4, 17.1, 17.2, 17.0, 17.3]
    
    # Each measurement has some uncertainty from the linear fit
    # (these would come from np.polyfit covariance matrix in real code)
    measurement_uncertainties = [0.020, 0.025, 0.018, 0.022, 0.019, 
                                0.023, 0.021, 0.020, 0.024, 0.021]
    
    n_grooves = len(true_angles)
    
    # OLD METHOD (what you had before):
    print("\n" + "="*70)
    print("OLD METHOD (Before Priority 1)")
    print("="*70)
    old_mean = np.mean(true_angles)
    old_std = np.std(true_angles, ddof=1)
    old_sem = old_std / np.sqrt(n_grooves)
    
    print(f"\nMean angle: {old_mean:.3f}°")
    print(f"Std deviation: {old_std:.3f}° (groove-to-groove variation)")
    print(f"SEM (naive): {old_sem:.3f}°")
    print(f"95% CI (naive): ±{1.96*old_sem:.3f}°")
    print("\nPROBLEM: This ignores measurement uncertainty!")
    print("         SEM is underestimated because it only accounts for")
    print("         physical variation, not measurement precision.")
    
    # NEW METHOD (with Priority 1 improvements):
    print("\n" + "="*70)
    print("NEW METHOD (After Priority 1)")
    print("="*70)
    
    # Separate the two sources of uncertainty
    mean_angle = np.mean(true_angles)
    
    # 1. Measurement uncertainty (from fits)
    measurement_variance = np.mean([u**2 for u in measurement_uncertainties])
    mean_meas_unc = np.mean(measurement_uncertainties)
    
    # 2. Physical variation (groove-to-groove)
    physical_variance = np.var(true_angles, ddof=1)
    physical_std = np.sqrt(physical_variance)
    
    # 3. Combined
    total_variance = physical_variance + measurement_variance
    total_std = np.sqrt(total_variance)
    sem_correct = total_std / np.sqrt(n_grooves)
    
    print(f"\nMean angle: {mean_angle:.3f}°")
    print(f"\nUncertainty Decomposition:")
    print(f"  1. Measurement uncertainty (avg per groove): {mean_meas_unc:.3f}°")
    print(f"     └─ Variance contribution: {measurement_variance:.6f}°²")
    print(f"  2. Physical variation (groove-to-groove):   {physical_std:.3f}°")
    print(f"     └─ Variance contribution: {physical_variance:.6f}°²")
    print(f"  3. Total combined uncertainty:              {total_std:.3f}°")
    print(f"     └─ Total variance: {total_variance:.6f}°²")
    print(f"\nStandard error of mean (correct): {sem_correct:.3f}°")
    print(f"95% confidence interval (correct): ±{1.96*sem_correct:.3f}°")
    print(f"\nRESULT: {mean_angle:.3f}° ± {1.96*sem_correct:.3f}° (95% CI)")
    
    # Comparison
    print("\n" + "="*70)
    print("COMPARISON")
    print("="*70)
    difference = sem_correct - old_sem
    percent_increase = (difference / old_sem) * 100
    
    print(f"\nOld SEM (ignoring measurement uncertainty): {old_sem:.3f}°")
    print(f"New SEM (including measurement uncertainty): {sem_correct:.3f}°")
    print(f"Difference: +{difference:.3f}° ({percent_increase:.1f}% larger)")
    print(f"\nOld 95% CI width: ±{1.96*old_sem:.3f}°")
    print(f"New 95% CI width: ±{1.96*sem_correct:.3f}°")
    
    print("\n" + "="*70)
    print("KEY INSIGHT")
    print("="*70)
    print("\nIn this example:")
    print(f"  • Each groove can be measured to ±{mean_meas_unc:.3f}° precision")
    print(f"  • Grooves vary by ±{physical_std:.3f}° due to real physical differences")
    print(f"  • The physical variation dominates (physical >> measurement)")
    print(f"  • But ignoring measurement uncertainty underestimates total uncertainty by {percent_increase:.0f}%")
    print("\nWith better AFM data (sharper features, less noise):")
    print("  • Measurement uncertainty would be even smaller")
    print("  • Difference between old and new methods would be smaller")
    print("\nWith worse AFM data (noisy, rounded features):")
    print("  • Measurement uncertainty would be larger")
    print("  • Difference between old and new methods would be larger")
    
    # Show what R² tells you about fit uncertainty
    print("\n" + "="*70)
    print("RELATIONSHIP TO R² VALUES")
    print("="*70)
    print("\nFit quality (R²) relates to measurement uncertainty:")
    print("  R² > 0.999  →  Very precise fit, uncertainty ~0.01-0.02°")
    print("  R² = 0.995  →  Good fit, uncertainty ~0.02-0.03°")
    print("  R² = 0.990  →  Acceptable fit, uncertainty ~0.03-0.05°")
    print("  R² < 0.980  →  Poor fit, uncertainty >0.05°, consider excluding")
    print("\nThe covariance matrix from np.polyfit quantifies this automatically!")


def show_sample_output():
    """Show what the new output looks like"""
    
    print("\n\n" + "="*70)
    print("SAMPLE OUTPUT FROM NEW CODE")
    print("="*70)
    
    sample_output = """
============================================================
RESULTS FOR sample_grating.txt
============================================================
Analysis side: negative_slope
Grooves analyzed: 10

Per-groove statistics (groove-to-groove variation):
  Mean blaze angle: 17.15 deg ± 0.18 deg (physical variation)
  Min/Max blaze angle: 16.90 deg / 17.40 deg
  Mean slope (dy/dx): 0.3087

Uncertainty analysis:
  Average measurement uncertainty per groove: 0.023 deg
  Physical variation (groove-to-groove): 0.180 deg
  Total uncertainty (combined): 0.181 deg
  Standard error of mean: 0.057 deg
  95% confidence interval on mean: ±0.112 deg

Within-facet statistics (camber/curvature):
  Local angle std: 0.15 deg
  Local angle range: 0.52 deg
  Mean within-facet variation: 0.15 deg

Groove geometry:
  Measured groove spacing: 833.33 nm ± 12.45 nm
  Mean groove depth: 245.67 nm
  Mean blaze facet width: 567.89 nm
"""
    
    print(sample_output)
    print("\n" + "="*70)
    print("WHAT THIS TELLS YOU")
    print("="*70)
    print("""
1. The mean blaze angle is 17.15°
   
2. Grooves vary by ±0.18° from each other (real physical differences)

3. Each individual measurement has ±0.023° uncertainty
   (much smaller than the physical variation - good!)

4. The true mean angle is between 17.04° and 17.26° (95% confidence)

5. Within each facet, the angle varies by ±0.15° (surface curvature)

6. All uncertainties are now properly quantified and separated!
""")


if __name__ == "__main__":
    demonstrate_uncertainty_decomposition()
    show_sample_output()
    
    print("\n" + "="*70)
    print("NEXT STEPS")
    print("="*70)
    print("""
1. Replace your analysis.py and analyzer.py with the updated versions
2. Run your normal analysis workflow
3. Check the new "Uncertainty analysis" section in output
4. Verify that measurement uncertainty < physical variation (usually 5-10x)
5. Use the 95% CI for reporting results in papers

Then we can move on to Priority 2: Visualization improvements!
""")
