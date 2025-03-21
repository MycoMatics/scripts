#!/bin/bash
# Script: check_inode_usage_with_log.sh
# Description:
#   This script recursively calculates inode usage for a specified directory
#   and its subdirectories. It logs each directory's inode count to a log file,
#   prints the grand total for the target directory, and at the end, displays
#   the 10 directories with the highest inode usage.
#
# Usage: ./check_inode_usage_with_log.sh <directory>

# Check if a target directory is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

TARGET_DIR="$1"

# Verify that the provided argument is a directory
if [ ! -d "$TARGET_DIR" ]; then
    echo "Error: Directory '$TARGET_DIR' does not exist."
    exit 1
fi

# Define the log file (timestamped for uniqueness)
LOG_FILE="inode_usage_$(date +%Y%m%d_%H%M%S).log"

# Write header to log file
{
    echo "Inode Usage Report for '$TARGET_DIR'"
    echo "Generated on $(date)"
    echo "---------------------------------------------------------"
} | tee "$LOG_FILE"

# Temporary file to store directory inode counts for sorting later
TMP_FILE=$(mktemp)

# Iterate over each directory recursively and compute inode usage.
# The -print0 option and IFS handling ensure correct processing of names with spaces.
find "$TARGET_DIR" -type d -print0 | while IFS= read -r -d '' dir; do
    # Count inodes in the current directory tree (files and directories)
    inode_count=$(find "$dir" -mindepth 0 | wc -l)
    output="Directory: $dir -> Inode count: $inode_count"
    
    # Print the result and append it to the log file
    echo "$output" | tee -a "$LOG_FILE"
    
    # Save the inode count and directory name to the temporary file for later sorting
    echo "$inode_count $dir" >> "$TMP_FILE"
done

# Calculate the grand total inode count for the entire target directory tree.
grand_total=$(find "$TARGET_DIR" -mindepth 0 | wc -l)
{
    echo "---------------------------------------------------------"
    echo "Grand Total for '$TARGET_DIR': $grand_total inodes"
    echo "---------------------------------------------------------"
} | tee -a "$LOG_FILE"

# Display the top 10 directories with the highest inode usage.
echo "Top 10 directories with highest inode usage:" | tee -a "$LOG_FILE"
echo "---------------------------------------------------------" | tee -a "$LOG_FILE"
# Sort the temporary file in descending numerical order by inode count,
# then print the top 10 entries.
sort -nr "$TMP_FILE" | head -n 10 | while read -r count dir; do
    echo "Directory: $dir -> Inode count: $count" | tee -a "$LOG_FILE"
done

# Remove the temporary file
rm "$TMP_FILE"

echo "Report saved to $LOG_FILE"
