#!/bin/bash

fullname=$(basename $1)
filename1=${fullname%_*.fq}
filename=${fullname%%_*}
echo ${fullname}
echo ${filename1}
echo ${filename}

------------EXAMPLE------------------------


$./name.sh sample1_cleaned_1.fq
sample1_cleaned_1.fq
sample1_cleaned
sample1
-------------------------------------------
Substring removal

${PARAMETER#PATTERN}

${PARAMETER##PATTERN}

${PARAMETER%PATTERN}

${PARAMETER%%PATTERN}

This one can expand only a part of a parameter's value, given a pattern to describe what to remove from the string. The pattern is interpreted just like a pattern to describe a filename to match (globbing). See Pattern matching for more.

Example string (just a quote from a big man):

MYSTRING="Be liberal in what you accept, and conservative in what you send"

From the beginning

${PARAMETER#PATTERN} and ${PARAMETER##PATTERN}

This form is to remove the described pattern trying to match it from the beginning of the string. The operator "#" will try to remove the shortest text matching the pattern, while "##" tries to do it with the longest text matching. Look at the following examples to get the idea (matched text marked striked, remember it will be removed!):
Syntax	Result
${MYSTRING#*in}	Be liberal in what you accept, and conservative in what you send
${MYSTRING##*in}	Be liberal in what you accept, and conservative in what you send
From the end

${PARAMETER%PATTERN} and ${PARAMETER%%PATTERN}

In the second form everything will be the same, except that Bash now tries to match the pattern from the end of the string:
Syntax 	Result
${MYSTRING%in*} 	Be liberal in what you accept, and conservative in what you send
${MYSTRING%%in*} 	Be liberal in what you accept, and conservative in what you send
