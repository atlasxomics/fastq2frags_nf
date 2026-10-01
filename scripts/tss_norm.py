import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
import logging

def normalize_tss_profile(tss_profile, window_size=1000, flank_size=100):
    """
    Normalize TSS profile using background normalization method.
    
    Parameters:
    -----------
    tss_profile : array-like
        Raw TSS profile signal (typically 2001 bp: -1000 to +1000 around TSS)
    window_size : int
        Size of window around TSS in bp (default 1000)
    flank_size : int
        Size of flanking regions for background calculation (default 100)
        
    Returns:
    --------
    normalized_profile : numpy.array
        Background-normalized TSS profile (fold enrichment)
    background_mean : float
        Mean signal in background regions
    """
    
    tss_profile = np.array(tss_profile)
    
    # Calculate indices for flanking regions
    # Assuming profile is centered at TSS (middle position)
    profile_length = len(tss_profile)
    center = profile_length // 2
    
    # Left flank: start to start+flank_size
    left_start = 0
    left_end = flank_size
    
    # Right flank: end-flank_size to end
    right_start = profile_length - flank_size
    right_end = profile_length
    
    # Extract flanking regions
    left_flank = tss_profile[left_start:left_end]
    right_flank = tss_profile[right_start:right_end]
    
    # Calculate background mean from flanking regions
    background_signal = np.concatenate([left_flank, right_flank])
    background_mean = np.mean(background_signal)
    
    # Avoid division by zero
    if background_mean == 0:
        background_mean = 1e-6
    
    # Normalize profile by background
    normalized_profile = tss_profile / background_mean
    
    return normalized_profile, background_mean

def calculate_tss_enrichment_score(normalized_profile, center_window=200):
    """
    Calculate TSS enrichment score from normalized profile.
    
    Parameters:
    -----------
    normalized_profile : array-like
        Background-normalized TSS profile
    center_window : int
        Window around TSS center to calculate max enrichment (default 200bp)
        
    Returns:
    --------
    tss_score : float
        Maximum fold enrichment within center window
    """
    
    profile_length = len(normalized_profile)
    center = profile_length // 2
    
    # Define window around TSS center
    window_start = max(0, center - center_window // 2)
    window_end = min(profile_length, center + center_window // 2)
    
    # Get maximum enrichment in center window
    center_region = normalized_profile[window_start:window_end]
    tss_score = np.max(center_region)
    
    return tss_score

def smooth_profile(profile, sigma=1):
    """Apply gaussian smoothing to profile."""
    return ndimage.gaussian_filter1d(profile, sigma=sigma)

# Example usage and demonstration
def demo_tss_normalization(tss_profile=None, window_size=1000, flank_size=100):
    """
    Demonstrate TSS profile normalization with provided data or synthetic data.
    
    Parameters:
    -----------
    tss_profile : array-like, optional
        TSS profile data to normalize. If None, generates synthetic data.
    window_size : int
        Size of window around TSS in bp (default 1000)
    flank_size : int
        Size of flanking regions for background calculation (default 100)
        
    Returns:
    --------
    tuple : (raw_profile, normalized_profile, positions)
        Raw profile, normalized profile, and position array
    """
    logging.info("="*60)
    logging.info("TSS NORMALIZATION PLOT")
    logging.info("="*60)
    if tss_profile is None:
        # Create synthetic TSS profile data
        logging.warning("No TSS profile provided. Generating synthetic data...")
        np.random.seed(42)
        
        # Base signal level
        base_level = 10
        
        # Create position array based on window size
        positions = np.arange(-window_size, window_size + 1)
        
        # Simulate TSS enrichment peak (gaussian-like)
        tss_peak = 50 * np.exp(-(positions**2) / (2 * 150**2))
        
        # Add some noise
        noise = np.random.normal(0, 2, len(positions))
        
        # Combine to create raw profile
        raw_profile = base_level + tss_peak + noise
        
        # Ensure no negative values
        raw_profile = np.maximum(raw_profile, 0.1)
        
        logging.info(f"Generated synthetic TSS profile with {len(raw_profile)} positions")
        
    else:
        # Use provided TSS profile
        raw_profile = np.array(tss_profile)
        profile_length = len(raw_profile)
        
        # Infer window size from profile length
        # Assume profile is symmetric around TSS
        inferred_window = (profile_length - 1) // 2
        positions = np.arange(-inferred_window, inferred_window + 1)
        
        logging.info(f"Using provided TSS profile with {len(raw_profile)} positions")
        logging.info(f"Inferred window size: ±{inferred_window} bp")

        # need to allow user to set window size
        if window_size is not None:
            center = (profile_length-1) // 2
            windowed_start = center - window_size # 2000-1000 = 1000 -1
            windowed_end = center + window_size +1 # 2000+1000 = 3000
            print(profile_length, center, windowed_start, windowed_end)
            windowed_raw_profile = raw_profile[windowed_start:windowed_end]
            windowed_profile_length = len(windowed_raw_profile)
            positions = np.arange(-window_size, window_size + 1)
            logging.info(f"\tWindowed TSS profile with {windowed_profile_length} positions")
            logging.info(f"\tRequested window size: ±{window_size} bp")
            raw_profile = windowed_raw_profile
        
        # Adjust window_size and flank_size if they're larger than profile
        if window_size > inferred_window:
            window_size = inferred_window
            logging.info(f"Adjusted window size to {window_size} bp (profile limit)")
            
        if flank_size * 2 >= profile_length:
            flank_size = max(1, profile_length // 10)  # Use 10% of profile length
            logging.info(f"Adjusted flank size to {flank_size} bp (profile limit)")
    
    # Normalize the profile
    norm_profile, bg_mean = normalize_tss_profile(raw_profile, 
                                                  window_size=window_size, 
                                                  flank_size=flank_size)
    
    # Calculate TSS enrichment score
    tss_score = calculate_tss_enrichment_score(norm_profile)
    
    # Smooth for visualization
    smooth_norm = smooth_profile(norm_profile, sigma=3)

    # Calculate TSS enrichment score
    tss_score = calculate_tss_enrichment_score(smooth_norm)
    
    # Create visualization
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
    
    # Calculate background region positions for visualization
    profile_length = len(raw_profile)
    left_flank_pos = positions[:flank_size]
    right_flank_pos = positions[-flank_size:]
    
    # Plot raw profile
    ax1.plot(positions, raw_profile, 'b-', alpha=0.7, label='Raw Profile')
    ax1.axhline(y=bg_mean, color='r', linestyle='--', 
                label=f'Background Mean = {bg_mean:.2f}')
    
    # Highlight background regions
    if len(left_flank_pos) > 0:
        ax1.axvspan(left_flank_pos[0], left_flank_pos[-1], 
                   alpha=0.2, color='red', label='Background Regions')
    if len(right_flank_pos) > 0:
        ax1.axvspan(right_flank_pos[0], right_flank_pos[-1], 
                   alpha=0.2, color='red')
    
    ax1.set_xlabel('Distance from TSS (bp)')
    ax1.set_ylabel('Raw Signal')
    ax1.set_title('Raw TSS Profile')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot normalized profile
    ax2.plot(positions, norm_profile, 'g-', alpha=0.7, label='Normalized')
    ax2.plot(positions, smooth_norm, 'k-', linewidth=2, label='Smoothed')
    ax2.axhline(y=1.0, color='r', linestyle='--', alpha=0.7, 
                label='Background Level (1.0)')
    
    # Highlight TSS center region
    # print(window_size // 2)
    center_window = min(200, window_size // 2)  # Adjust center window to profile size
    ax2.axvspan(-center_window // 2, center_window // 2, alpha=0.2, color='yellow', 
                label='TSS Center Region')
    
    ax2.set_xlabel('Distance from TSS (bp)')
    ax2.set_ylabel('Fold Enrichment over Background')
    ax2.set_title(f'Normalized TSS Profile (TSS Score = {tss_score:.2f})')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.close()
    #plt.show()
    
    # Print results
    logging.info("TSS Normalization Results:")
    # print("=" * 30)
    logging.info(f"Profile length:\t{len(raw_profile)} positions")
    logging.info(f"Window size:\t±{len(positions)//2} bp")
    logging.info(f"Flank size:\t{flank_size} bp")
    logging.info(f"Background mean signal:\t{bg_mean:.3f}")
    logging.info(f"TSS enrichment score:\t{tss_score:.3f}")
    logging.info(f"Max fold enrichment:\t{np.max(norm_profile):.3f}")
    logging.info(f"Min fold enrichment:\t{np.min(norm_profile):.3f}")
    
    return fig, raw_profile, norm_profile, positions

# Function to work with real SnapATAC2 data
def normalize_snapatac2_tsse_profile(adata, cell_id=None):
    """
    Normalize TSS profile from SnapATAC2 AnnData object.
    
    Parameters:
    -----------
    adata : anndata.AnnData
        SnapATAC2 AnnData object with TSSE calculated
    cell_id : str, optional
        Specific cell ID to normalize (if None, uses aggregate)
        
    Returns:
    --------
    normalized_profile : numpy.array
        Normalized TSS profile
    """
    
    # Extract TSS profile from adata
    # This assumes the profile is stored in adata.uns or similar
    # Adjust the key based on your SnapATAC2 version
    
    if 'tsse_profile' in adata.uns:
        raw_profile = adata.uns['tsse_profile']
    elif hasattr(adata, 'obsm') and 'tsse_profile' in adata.obsm:
        if cell_id:
            cell_idx = list(adata.obs.index).index(cell_id)
            raw_profile = adata.obsm['tsse_profile'][cell_idx]
        else:
            # Use mean across all cells
            raw_profile = np.mean(adata.obsm['tsse_profile'], axis=0)
    else:
        raise ValueError("TSS profile not found. Make sure TSSE has been calculated.")
    
    # Normalize using background method
    normalized_profile, bg_mean = normalize_tss_profile(raw_profile)
    
    return normalized_profile, bg_mean

if __name__ == "__main__":
    # Run demonstration with synthetic data
    print("Demonstrating TSS Profile Background Normalization")
    print("=" * 50)
    demo_tss_normalization()
    
    print("\n" + "=" * 50)
    print("Example with custom TSS profile:")
    
    # Example with a shorter custom profile
    custom_profile = [5, 6, 8, 12, 25, 45, 60, 55, 30, 15, 10, 8, 6, 5, 4]
    tfig, _, _, _ = demo_tss_normalization(custom_profile, flank_size=200)
    