#!/usr/bin/env python3
"""
Fragment Analysis Report Generator

This script generates an HTML report for fragment analysis results,
embedding images and extracting metrics from log files.

Usage:
    python generate_report.py <directory_path>

Requirements:
    - PIL (Pillow) for image processing
    - Standard library modules: os, sys, re, base64, datetime
"""

import os
import sys
import re
import base64
import datetime
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional

try:
    from PIL import Image
    PIL_AVAILABLE = True
except ImportError:
    PIL_AVAILABLE = False
    print("Warning: PIL (Pillow) not found. Image optimization will be skipped.")


def encode_image_to_base64(image_path: str) -> Optional[str]:
    """
    Convert an image file to a base64 data URL string.
    
    Args:
        image_path: Path to the image file
        
    Returns:
        Base64 data URL string or None if file doesn't exist
    """
    if not os.path.exists(image_path):
        print(f"Warning: Image file not found: {image_path}")
        return None
        
    try:
        # Determine MIME type based on file extension
        ext = Path(image_path).suffix.lower()
        mime_types = {
            '.png': 'image/png',
            '.jpg': 'image/jpeg',
            '.jpeg': 'image/jpeg',
            '.gif': 'image/gif',
            '.svg': 'image/svg+xml'
        }
        mime_type = mime_types.get(ext, 'image/png')
        
        # Read and encode the image
        with open(image_path, 'rb') as img_file:
            img_data = img_file.read()
            
        # Optionally compress PNG images if PIL is available
        if PIL_AVAILABLE and ext == '.png':
            try:
                img = Image.open(image_path)
                # Convert to RGB if necessary (removes alpha channel)
                if img.mode in ('RGBA', 'LA'):
                    background = Image.new('RGB', img.size, (255, 255, 255))
                    if img.mode == 'RGBA':
                        background.paste(img, mask=img.split()[-1])
                    else:
                        background.paste(img, mask=img.split()[-1])
                    img = background
                
                # Save optimized version to bytes
                from io import BytesIO
                buffer = BytesIO()
                img.save(buffer, format='PNG', optimize=True)
                img_data = buffer.getvalue()
                print(f"Optimized image: {image_path}")
            except Exception as e:
                print(f"Could not optimize {image_path}: {e}")
        
        encoded = base64.b64encode(img_data).decode('utf-8')
        return f"data:{mime_type};base64,{encoded}"
        
    except Exception as e:
        print(f"Error encoding image {image_path}: {e}")
        return None


def parse_chromap_log_stats(log_file: str, req_string: str) -> float:
    """
    Extract specific metrics from chromap log file.
    
    Args:
        log_file: Path to the log file
        req_string: String to search for in the log file
        
    Returns:
        The value associated with the string
    """
    # Read lines and try to find the requested string. Return a parsed value
    # which may be int/float or a string. This function is tolerant to
    # different separators and will fall back to returning the raw string
    # if parsing fails.
    if not os.path.exists(log_file):
        return ""

    with open(log_file, "r", encoding="utf-8", errors="ignore") as file:
        for line in file:
            if req_string in line:
                # Try to extract the text after the first ':' if present,
                # otherwise take the portion after the matched req_string.
                try:
                    if ':' in line:
                        raw = line.split(':', 1)[1].strip()
                    else:
                        raw = line.split(req_string, 1)[1].strip()
                except Exception:
                    raw = line.strip()

                # Strip common trailing punctuation
                raw = raw.rstrip('\n').strip()
                if raw.endswith('.') or raw.endswith(','):
                    raw = raw[:-1].strip()

                # No structured parsing: this log only uses plain numbers or strings.
                # Proceed to numeric parsing below.
                try:
                    if re.match(r'^[-+]?[0-9]+\.[0-9]+$', raw):
                        return float(raw)
                    if re.match(r'^[-+]?[0-9]+$', raw):
                        return int(raw)
                except Exception:
                    pass

                # Strip surrounding quotes if present and return raw string
                if (raw.startswith('"') and raw.endswith('"')) or (raw.startswith("'") and raw.endswith("'")):
                    return raw[1:-1]
                return raw

    # If nothing found, return empty string
    return ""

def parse_chromap_metrics(log_path: str, singlecell_df: pd.DataFrame = None) -> List[Dict]:
    """
    Parse chromap-specific metrics from the log file.
    
    Args:
        log_path: Path to the log file
        singlecell_df: Optional DataFrame with single cell data
        
    Returns:
        List of metrics
    """
    metrics = []

    def _to_number(x):
        """Attempt to coerce x to a float. Return None if not possible."""
        if x is None:
            return None
        if isinstance(x, (int, float)):
            return float(x)
        if isinstance(x, str):
            # Remove commas and percent signs
            s = x.replace(',', '').replace('%', '').strip()
            # Remove trailing text like 'read pairs' if present by taking first token
            # but only if it's clearly numeric-looking
            try:
                return float(s)
            except Exception:
                # Try to extract leading number
                m = re.match(r'^[\+\-]?[0-9]*\.?[0-9]+', s)
                if m:
                    try:
                        return float(m.group(0))
                    except Exception:
                        return None
                return None
        return None

    try:
        # Raw values from chromap log
        uni_map = parse_chromap_log_stats(log_path, "Number of uni-mappings")
        mappings = parse_chromap_log_stats(log_path, "Number of mappings")
        mapped_reads = parse_chromap_log_stats(log_path, "Number of mapped reads")
        total_reads = parse_chromap_log_stats(log_path, "Number of reads")
        corrected_bcs = parse_chromap_log_stats(log_path, "Number of corrected barcodes")
        barcodes_in_whitelist = parse_chromap_log_stats(log_path, "Number of barcodes in whitelist")

        # Add raw counts if available
        for key, raw in (
            ("Number of uni-mappings", uni_map),
            ("Number of mappings", mappings),
            ("Number of mapped reads", mapped_reads),
            ("Number of reads", total_reads),
            ("Number of corrected barcodes", corrected_bcs),
            ("Number of barcodes in whitelist", barcodes_in_whitelist),
        ):
            if raw not in (None, ""):
                metrics.append({
                    'name': key,
                    'value': str(raw),
                    'description': 'Raw value extracted from chromap log'
                })

        # Compute fractions only when numeric
        uni_num = _to_number(uni_map)
        map_num = _to_number(mappings)
        if uni_num is not None and map_num:
            try:
                frac_uniq = uni_num / map_num if map_num != 0 else None
                if frac_uniq is not None:
                    metrics.append({
                        'name': 'Fraction Uniquely Aligned Reads',
                        'value': f"{frac_uniq:.4f}",
                        'description': 'Proportion of reads that aligned uniquely to the genome'
                    })
            except Exception:
                pass

        mapped_num = _to_number(mapped_reads)
        total_num = _to_number(total_reads)
        if total_num:
            try:
                metrics.append({
                    'name': 'Number of input read pairs',
                    'value': f"{float(total_num/2):.1f}",
                    'description': 'Total number input read pairs (assuming paired-end reads)'
                })
            except Exception:
                pass

        if mapped_num is not None and total_num:
            try:
                frac_unaligned = 1 - (mapped_num / total_num) if total_num != 0 else None
                if frac_unaligned is not None:
                    metrics.append({
                        'name': 'Fraction Unaligned Reads',
                        'value': f"{frac_unaligned:.4f}",
                        'description': 'Proportion of reads that failed to align'
                    })
            except Exception:
                pass

        corr_num = _to_number(corrected_bcs)
        whitelist_num = _to_number(barcodes_in_whitelist)
        if corr_num is not None and whitelist_num:
            try:
                frac_valid_bc = 1 - (corr_num / whitelist_num) if whitelist_num != 0 else None
                if frac_valid_bc is not None:
                    metrics.append({
                        'name': 'Fraction Reads with Valid Barcode',
                        'value': f"{frac_valid_bc:.4f}",
                        'description': 'Proportion of reads with valid cell barcodes'
                    })
            except Exception:
                pass
        
        # If singlecell dataframe provided, report total read pairs from it
        if singlecell_df is not None and 'total' in singlecell_df.columns:
            try:
                total_from_df = int(singlecell_df['total'].sum())
                metrics.append({
                    'name': 'Chromap Input Read Pairs',
                    'value': f"{total_from_df}",
                    'description': 'Total number of read pairs processed (from singlecell.csv)'
                })
            except Exception:
                pass

        # Additional string-like entries
        bcListFile = parse_chromap_log_stats(log_path, "Cell barcode whitelist file")
        if bcListFile not in (None, ""):
            metrics.append({
                'name': 'Cell Barcode Whitelist File',
                'value': str(bcListFile),
                'description': 'Path to the cell barcode whitelist file used by Chromap'
            })

        # Sample run id from run_id.txt if present
        try:
            log_parent_path = Path(log_path).parent.resolve()
            run_id_path = log_parent_path.joinpath("run_id.txt")
            sampleRun_id = open(run_id_path, "r").read().strip()
            if sampleRun_id:
                metrics.append({
                    'name': 'Sample Run ID',
                    'value': sampleRun_id,
                    'description': 'Unique identifier for the sample run'
                })
        except Exception:
            print(f"No run_id.txt file found: {run_id_path}")
            pass

        # Workflow version if present
        try:
            version_path = Path("/root/version").resolve()
            version_id = open(version_path, "r").read().strip()
            if version_id:
                metrics.append({
                    'name': 'Atlas Pipeline Version',
                    'value': version_id,
                    'description': 'Version running'
                })
        except Exception:
            print(f"No version file found: {version_path}")
            pass

    except Exception as e:
        print(f"Error parsing chromap metrics: {e}")

    return metrics

def parse_log_file(log_path: str) -> Tuple[List[Dict], str]:
    """
    Parse the log file to extract metrics and preview content.
    
    Args:
        log_path: Path to the log file
        
    Returns:
        Tuple of (metrics_list, preview_content)
    """
    metrics = []
    preview_content = ""
    
    if not os.path.exists(log_path):
        print(f"Warning: Log file not found: {log_path}")
        return metrics, preview_content
    
    try:
        with open(log_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
            lines = content.split('\n')
        
        # Generate preview (first 20 lines)
        preview_lines = lines[:20]
        preview_content = '\n'.join(preview_lines)
        
        # Define patterns for fragment analysis log format (based on the provided example)
        patterns = [
            # Basic metrics
            (r'Total cells analyzed:\s+(\d+)', 'Total Barcodes Analyzed', 'Total number of barcodes passing filters'),
            (r'mean nfrag:\s+([\d.]+)', 'Mean Fragments per pixel', 'Average number of fragments per pixel'),
            (r'median nfrag:\s+([\d.]+)', 'Median Fragments per pixel', 'Median number of fragments per pixel'),
            
            # Duplication metrics
            (r'mean dup:\s+([\d.e-]+)', 'Mean Duplication Rate', 'Average duplication rate per pixel'),
            (r'median dup:\s+([\d.e-]+)', 'Median Duplication Rate', 'Median duplication rate per pixel'),
            (r'frag_wt_dup:\s+([\d.e-]+)', 'Fragment-weighted Duplication', 'Fragment-weighted duplication rate'),
            
            # Mitochondrial metrics
            (r'mito_mean:\s+([\d.]+)', 'Mean Mitochondrial Fraction', 'Average mitochondrial fragment fraction'),
            (r'mito_median:\s+([\d.]+)', 'Median Mitochondrial Fraction', 'Median mitochondrial fragment fraction'),
            
            # TSS metrics
            (r'library_tsse:\s+([\d.]+)', 'Library TSS Enrichment', 'TSS enrichment score (SNAPATAC2)'),
            (r'frac_overlap_TSS:\s+([\d.]+)', 'TSS Overlap Fraction', 'Fraction of fragments overlapping TSS'),
            (r'TSS enrichment score:\s+([\d.]+)', 'TSS Enrichment Score', 'TSS enrichment score (plot norm)'),
            (r'Max fold enrichment:\s+([\d.]+)', 'Max Fold Enrichment', 'Maximum fold enrichment at TSS'),
            (r'Min fold enrichment:\s+([\d.]+)', 'Min Fold Enrichment', 'Minimum fold enrichment at TSS'),
            (r'Background mean signal:\s+([\d.]+)', 'Background Mean Signal', 'Mean background signal level'),

            # FRiP metrics (500 bp merge_peaks consensus, blacklist-filtered, MACS3 q<=0.05;
            # peak-calling barcode set depends on --input-fragment-floor, see fragments_log.txt)
            (r'mean frip:\s+([\d.e-]+)', 'Mean FRiP', 'Average Fraction Reads in Peaks (500 bp blacklist-filtered consensus)'),
            (r'median frip:\s+([\d.e-]+)', 'Median FRiP', 'Median Fraction Reads in Peaks (500 bp blacklist-filtered consensus)'),
            (r'frag_wt_frip:\s+([\d.e-]+)', 'Fragment-weighted FRiP', 'Fragment-weighted Fraction Reads in Peaks (500 bp blacklist-filtered consensus)'),
            (r'frip_calibrated:\s+(NA|[\d.e-]+)', 'FRiP (calibrated)', 'FRiP calibrated to the pycisTopic reference (mode-selected fit; constants logged as frip_calib_params; NA when withheld by the coverage-collapse gate, see FRiP Reliability Flags)'),
            (r'frip_reliable:\s+(\w+)', 'FRiP Calibration Reliable', 'Whether per-barcode coverage supports the calibrated FRiP; False => the calibrated value is withheld (NA)'),
            (r'frip_flags:\s+([\w;]+)', 'FRiP Reliability Flags', 'Coverage-collapse reasons the calibration was withheld (none when reliable): calc_floor_fallback; low_frac_pass; low_median_nfrag'),
            (r'frac_pass:\s+([\d.e-]+)', 'Fraction Barcodes >= Calc Floor', 'Fraction of barcodes clearing the calc floor (min fragments); low values mean the FRiP calibration is out of scope'),
            (r'total_npeak:\s+([\d]+)', 'Total Peaks Called', 'Total number of called peaks'),
            
            # Profile metrics
            (r'Profile length:\s+(\d+)\s+positions', 'Profile Length', 'Number of positions in TSS profile'),
            (r'Window size:\s+Â±(\d+)\s+bp', 'Window Size', 'TSS analysis window size (±bp)'),
            (r'Flank size:\s+(\d+)\s+bp', 'Flank Size', 'Flanking region size for background calculation (bp)'),
            
            # Timing information
            (r'Fragment import:\s+([\d.]+)\s+(\w+)', 'Fragment Import Time', 'Time taken for fragment import'),
            (r'Plotting/metrics:\s+([\d.]+)\s+(\w+)', 'Analysis Time', 'Time taken for plotting and metrics'),
            (r'Bigwig export:\s+([\d.]+)\s+(\w+)', 'Bigwig Export Time', 'Time taken for bigwig export'),
            (r'FRiP calculation:\s+([\d.]+)\s+(\w+)', 'FRiP Calc Time', 'Time taken for FRiP calculation'),
            (r'Total runtime:\s+([\d.]+)\s+(\w+)', 'Total Runtime', 'Total analysis runtime'),
            
            # File information
            (r'Genome:\s+(\w+)', 'Genome Build', 'Reference genome used for analysis'),
        ]
        
        # Extract metrics using patterns
        for line in lines:
            for pattern, name, description in patterns:
                match = re.search(pattern, line, re.IGNORECASE)
                if match and not any(m['name'] == name for m in metrics):
                    value = match.group(1)
                    try:
                        units = match.group(2)
                        description = f"{description} ({units})"
                    except:
                        units = ""
                    metrics.append({
                        'name': name,
                        'value': value,
                        'description': description
                    })
        
        # Format scientific notation values for better readability
        for metric in metrics:
            try:
                value = float(metric['value'])
                if 'e-' in metric['value']:  # Scientific notation with negative exponent
                    if value < 0.001:
                        metric['value'] = f"{value:.2e}"
                    else:
                        metric['value'] = f"{value:.6f}"
                elif value > 1000:
                    metric['value'] = f"{value:,.0f}"
                elif value > 1:
                    metric['value'] = f"{value:.2f}"
                else:
                    metric['value'] = f"{value:.4f}"
            except (ValueError, TypeError):
                pass  # Keep original value if conversion fails
        
        # If no specific metrics found, try to extract general numeric values
        if not metrics:
            for line in lines:
                # Look for pattern: "label: number"
                match = re.search(r'^([a-zA-Z\s_-]+)[:\s]+([\d,]+(?:\.[\d]+)?)', line.strip())
                if match and len(metrics) < 15:
                    label = match.group(1).strip()
                    value = match.group(2)
                    if 3 <= len(label) <= 50 and label not in [m['name'] for m in metrics]:
                        metrics.append({
                            'name': label.title(),
                            'value': value,
                            'description': 'Extracted from log file'
                        })
        
        print(f"Extracted {len(metrics)} metrics from log file")
        
    except Exception as e:
        print(f"Error reading log file {log_path}: {e}")
    
    return metrics, preview_content


def generate_html_report(directory: str, output_path: str) -> None:
    """
    Generate the complete HTML report.
    
    Args:
        directory: Directory containing the fragment analysis files
        output_path: Path where the HTML report should be saved
    """

    def _format_number(value):
        """Format a number (or numeric string) with commas as thousands separators."""
        try:
            return f"{int(value):,}"
        except (ValueError, TypeError):
            return value
    
    # Define expected files
    # image_files entries may be either 3-tuples: (filename, title, alt)
    # or 4-tuples: (filename, title, alt, cols) where cols is 1 or 2
    # If cols==2 the image card will span two grid columns.
    image_files = [
        ('fragments_dup_hist.png', 'Duplication Histogram', 'Fragment Duplication Histogram', 1),
        ('fragments_frag_size.png', 'Fragment Size Distribution', 'Fragment Size Distribution', 1),
        ('fragments_frag_v_dup.png', 'Fragment vs Duplication', 'Fragment vs Duplication Analysis', 1),
        ('fragments_mito_hist.png', 'Mitochondrial Histogram', 'Mitochondrial Fragment Histogram', 1),
        ('fragments_norm_tss.png', 'Normalized TSS', 'Normalized TSS Analysis', 1),
        ('fragments_frip_hist.png', 'FRiP Distribution', 'Fraction of Reads in Peaks', 1),
        # Make lane_qc span two columns by default
        ('lane_qc.png', 'Lane QC Plot', 'Lane Quality Control Plot', 2),
    ]
    # frag analysis log
    log_file = 'fragments_log.txt'
    
    # Parse log files and metrics
    log_path = os.path.join(directory, log_file)
    metrics, log_preview = parse_log_file(log_path)
    
    # Add chromap metrics if log file exists
    # chromap log
    chromap_log_path = os.path.join(directory, 'chromap_log.txt')
    if os.path.exists(chromap_log_path):
        singlecell_df = None
        chromap_metrics = parse_chromap_metrics(chromap_log_path, singlecell_df)
        metrics.extend(chromap_metrics)
    
    # save metrics as csv
    metrics_df = pd.DataFrame(metrics)
    # remove any commas used as thousands separators
    try:
        metrics_df = metrics_df.applymap(lambda x: str(x).replace(',', '') if isinstance(x, str) else x)
    except AttributeError as e:
        print(f"Error in .applymap, using .map: {e}")
        metrics_df = metrics_df.map(lambda x: str(x).replace(',', '') if isinstance(x, str) else x)

    metrics_csv_path = os.path.join(directory, 'summary.csv')
    print(f"Saving metrics to CSV at: {metrics_csv_path}")
    metrics_df.T.to_csv(metrics_csv_path, index=False, header=False)   
    
    sample_run_id = next((item['value'] for item in chromap_metrics if item['name'] == 'Sample Run ID'), None)
    # Process images
    image_cards_html = []
    for entry in image_files:
        # Support both 3- and 4-element tuples
        if len(entry) == 3:
            filename, title, alt_text = entry
            cols = 1
        else:
            filename, title, alt_text, cols = entry

        image_path = os.path.join(directory, filename)
        image_data = encode_image_to_base64(image_path)

        if image_data:
            img_html = f'<img src="{image_data}" alt="{alt_text}" style="max-width: 100%; height: auto; border-radius: 4px; border: 1px solid #e9ecef;">'
        else:
            img_html = f'''<div style="padding:40px; background:#f8f9fa; border:2px dashed #dee2e6; border-radius:4px; text-align:center; color:#6c757d;">
                Image not found: {filename}
            </div>'''

        span_class = ' span-2' if int(cols) == 2 else ''

        image_cards_html.append(f'''
            <div class="image-card{span_class}">
                <h3>{title}</h3>
                {img_html}
            </div>
        ''')
    
    # Generate metrics table HTML
    if metrics:
        metrics_rows = '\n'.join([
            f'''<tr>
                <td><strong>{metric['name']}</strong></td>
                <td class="metric-value">{_format_number(metric['value'])}</td>
                <td>{metric['description']}</td>
            </tr>''' for metric in metrics
        ])
        metrics_placeholder = ''
    else:
        metrics_rows = '''<tr>
            <td colspan="3" style="text-align: center; color: #6c757d; font-style: italic;">
                No metrics could be automatically extracted from the log file. Please review the log file manually.
            </td>
        </tr>'''
        metrics_placeholder = '''<div class="placeholder-message">
            <strong>Note:</strong> No metrics were automatically extracted from the log file. The log file may use a different format than expected.
        </div>'''
    
    # Generate log preview HTML
    if log_preview.strip():
        log_preview_html = f'<pre>{log_preview}</pre>'
        log_placeholder = ''
    else:
        log_preview_html = '<p style="color: #6c757d; font-style: italic;">Log file not found or empty.</p>'
        log_placeholder = '''<div class="placeholder-message">
            Log file <code>fragments_log.txt</code> was not found in the specified directory.
        </div>'''
    
    # Generate current timestamp
    current_time = datetime.datetime.now().strftime('%B %d, %Y at %I:%M %p')
    
    # HTML template
    html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Fragment Analysis Report</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            margin: 0;
            padding: 20px;
            background-color: #f8f9fa;
            color: #333;
        }}
        
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            overflow: hidden;
        }}
        
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            text-align: center;
        }}
        
        .header h1 {{
            margin: 0;
            font-size: 2.5em;
            font-weight: 300;
        }}
        
        .header p {{
            margin: 10px 0 0 0;
            opacity: 0.9;
            font-size: 1.1em;
        }}
        
        .content {{
            padding: 30px;
        }}
        
        .section {{
            margin-bottom: 40px;
        }}
        
        .section h2 {{
            color: #667eea;
            border-bottom: 2px solid #667eea;
            padding-bottom: 10px;
            margin-bottom: 20px;
            font-size: 1.8em;
        }}
        
        .metrics-table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        
        .metrics-table thead {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
        }}
        
        .metrics-table th,
        .metrics-table td {{
            padding: 15px;
            text-align: left;
            border-bottom: 1px solid #e9ecef;
        }}
        
        .metrics-table th {{
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            font-size: 0.9em;
        }}
        
        .metrics-table tbody tr:hover {{
            background-color: #f8f9fa;
            transition: background-color 0.3s ease;
        }}
        
        .metrics-table tbody tr:last-child td {{
            border-bottom: none;
        }}
        
        .metric-value {{
            font-weight: 600;
            color: #495057;
        }}
        
        .image-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        
        .image-card {{
            background: white;
            border-radius: 8px;
            padding: 20px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            text-align: center;
        }}
        
        .image-card h3 {{
            margin: 0 0 15px 0;
            color: #495057;
            font-size: 1.2em;
        }}
        
        /* Allow some image cards to span two columns when requested */
        .image-card.span-2 {{
            grid-column: span 2;
        }}

        /* On narrow screens, force spanning cards to occupy a single column */
        @media (max-width: 800px) {{
            .image-card.span-2 {{
                grid-column: span 1;
            }}
        }}
        
        .log-preview {{
            background: #f8f9fa;
            border: 1px solid #dee2e6;
            border-radius: 8px;
            padding: 20px;
            margin: 20px 0;
            font-family: 'Courier New', monospace;
            font-size: 0.9em;
            max-height: 400px;
            overflow-y: auto;
        }}
        
        .placeholder-message {{
            background: #fff3cd;
            border: 1px solid #ffeaa7;
            border-radius: 8px;
            padding: 15px;
            margin: 20px 0;
            color: #856404;
        }}
        
        .file-info {{
            background: #d1ecf1;
            border: 1px solid #bee5eb;
            border-radius: 8px;
            padding: 15px;
            margin: 20px 0;
            color: #0c5460;
        }}
        
        .copy-button {{
            background: linear-gradient(135deg, #28a745 0%, #20c997 100%);
            color: white;
            border: none;
            padding: 10px 20px;
            border-radius: 6px;
            cursor: pointer;
            font-size: 14px;
            font-weight: 600;
            margin: 10px 0;
            transition: all 0.3s ease;
            display: inline-flex;
            align-items: center;
            gap: 8px;
        }}
        
        .copy-button:hover {{
            background: linear-gradient(135deg, #218838 0%, #1e9e8a 100%);
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(40, 167, 69, 0.3);
        }}
        
        .copy-button:active {{
            transform: translateY(0);
        }}
        
        .copy-button.copied {{
            background: linear-gradient(135deg, #17a2b8 0%, #6f42c1 100%);
        }}
        
        .copy-icon {{
            width: 16px;
            height: 16px;
        }}
        
        .table-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 20px;
        }}
        
        .table-header h2 {{
            margin: 0;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>Fragment Analysis Report</h1>
            <p>Generated on {current_time}</p>
        </div>
        
        <div class="content">
            <div class="file-info">
                <strong>Sample ID:</strong> {sample_run_id}
            </div>
            
            <!-- Summary Metrics Section -->
            <div class="section">
                <div class="table-header">
                    <h2>📊 Summary Metrics</h2>
                    <button class="copy-button" onclick="copyMetricsToClipboard()" id="copyBtn">
                        <svg class="copy-icon" viewBox="0 0 24 24" fill="currentColor">
                            <path d="M16 1H4c-1.1 0-2 .9-2 2v14h2V3h12V1zm3 4H8c-1.1 0-2 .9-2 2v14c0 1.1.9 2 2 2h11c1.1 0 2-.9 2-2V7c0-1.1-.9-2-2-2zm0 16H8V7h11v14z"/>
                        </svg>
                        Copy to Excel
                    </button>
                </div>
                {metrics_placeholder}
                
                <table class="metrics-table">
                    <thead>
                        <tr>
                            <th>Metric</th>
                            <th>Value</th>
                            <th>Description</th>
                        </tr>
                    </thead>
                    <tbody>
                        {metrics_rows}
                    </tbody>
                </table>
            </div>
            
            <!-- Log File Preview Section -->
            <div class="section">
                <h2>📝 Log File Preview</h2>
                {log_placeholder}
                <div class="log-preview">
                    {log_preview_html}
                </div>
            </div>
            
            <!-- Analysis Charts Section -->
            <div class="section">
                <h2>📈 Analysis Charts</h2>
                <div class="image-grid">
                    {''.join(image_cards_html)}
                </div>
            </div>
        </div>
    </div>
    
    <script>
        async function copyMetricsToClipboard() {{
            const table = document.querySelector('.metrics-table');
            const button = document.getElementById('copyBtn');
            
            if (!table) {{
                alert('No metrics table found to copy');
                return;
            }}
            
            let csvContent = '';
            
            // Get headers
            const headers = table.querySelectorAll('thead th');
            const headerRow = Array.from(headers).map(th => th.textContent.trim()).join('\\t');
            csvContent += headerRow + '\\n';
            
            // Get data rows
            const dataRows = table.querySelectorAll('tbody tr');
            dataRows.forEach(row => {{
                const cells = row.querySelectorAll('td');
                if (cells.length > 0) {{
                    const rowData = Array.from(cells).map(td => {{
                        // Clean up the cell content
                        let text = td.textContent.trim();
                        // Remove any HTML formatting but keep the text
                        text = text.replace(/\\s+/g, ' ');
                        // Escape quotes and handle commas for CSV
                        if (text.includes('\\t') || text.includes('\\n') || text.includes('"')) {{
                            text = '"' + text.replace(/"/g, '""') + '"';
                        }}
                        return text;
                    }}).join('\\t');
                    csvContent += rowData + '\\n';
                }}
            }});
            
            try {{
                // Use the modern Clipboard API
                await navigator.clipboard.writeText(csvContent);
                
                // Visual feedback
                const originalContent = button.innerHTML;
                button.classList.add('copied');
                button.innerHTML = `
                    <svg class="copy-icon" viewBox="0 0 24 24" fill="currentColor">
                        <path d="M9 16.17L4.83 12l-1.42 1.41L9 19 21 7l-1.41-1.41z"/>
                    </svg>
                    Copied!
                `;
                
                setTimeout(() => {{
                    button.classList.remove('copied');
                    button.innerHTML = originalContent;
                }}, 2000);
                
            }} catch (err) {{
                // Fallback for older browsers
                const textArea = document.createElement('textarea');
                textArea.value = csvContent;
                textArea.style.position = 'fixed';
                textArea.style.left = '-999999px';
                textArea.style.top = '-999999px';
                document.body.appendChild(textArea);
                textArea.focus();
                textArea.select();
                
                try {{
                    document.execCommand('copy');
                    
                    // Visual feedback
                    const originalContent = button.innerHTML;
                    button.classList.add('copied');
                    button.innerHTML = `
                        <svg class="copy-icon" viewBox="0 0 24 24" fill="currentColor">
                            <path d="M9 16.17L4.83 12l-1.42 1.41L9 19 21 7l-1.41-1.41z"/>
                        </svg>
                        Copied!
                    `;
                    
                    setTimeout(() => {{
                        button.classList.remove('copied');
                        button.innerHTML = originalContent;
                    }}, 2000);
                    
                }} catch (err2) {{
                    console.error('Failed to copy metrics:', err2);
                    alert('Failed to copy to clipboard. Please try selecting the table manually.');
                }} finally {{
                    document.body.removeChild(textArea);
                }}
            }}
        }}
    </script>
</html>'''
    
    # Write the HTML file
    try:
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        print(f"Report generated successfully: {output_path}")
    except Exception as e:
        print(f"Error writing HTML file: {e}")
        sys.exit(1)


def main():
    """Main function to handle command line arguments and generate report."""
    if len(sys.argv) != 2:
        print("Usage: python generate_report.py <directory_path>")
        print("\nExample:")
        print("  python generate_report.py /path/to/fragment/analysis/results/")
        sys.exit(1)
    
    directory = sys.argv[1]
    
    # Validate directory
    if not os.path.isdir(directory):
        print(f"Error: Directory not found: {directory}")
        sys.exit(1)
    
    # Set output path
    output_path = os.path.join(directory, "fragment_analysis_report.html")

    print("Generating fragment analysis report...")
    print(f"Input directory: {directory}")
    print(f"Output file: {output_path}")
    print()
    
    # Generate the report
    generate_html_report(directory, output_path)
    
    # Show summary
    print()
    print("Report generation complete!")
    print("Open the following file in your browser to view the report:")
    print(f"  {output_path}")


if __name__ == "__main__":
    main()
