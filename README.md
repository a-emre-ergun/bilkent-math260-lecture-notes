# Bilkent MATH 260 Lecture Notes

Lecture notes for MATH 260 at Bilkent University, covering the final exam material of the Spring 2026 semester. I will add the notes for the earlier topics later.

I did not miss any lecture and wrote down everything on the board, so nothing from these lectures is missing. I had already turned my handwritten notes into a PDF for my convenience and wanted to share them. Claude Opus 5.5 was used to convert the notes into LaTeX and Markdown versions, and the mistakes in my handwritten notes were also fixed during the conversion. If you still find a mistake, please open an issue.

Keep in mind that these notes are from Spring 2026, so the topics may change in later semesters. I highly recommend studying from the textbook and solving its exercises, and using these notes as a quick reminder, since the textbook gives a much better understanding. Still, if you are short on time, the notes are a good way to learn the topics.

## Contents

### Final exam material

The chapters are in the order they were taught in the Spring 2026 semester. Chapter numbers follow the textbook, *Modern Mathematical Statistics with Applications* by Devore and Berk.

1. [Chapter 9: Tests of Hypotheses Based on a Single Sample](notes/ch09-hypothesis-tests-single-sample.md)
2. [Chapter 10: Inferences Based on Two Samples](notes/ch10-inferences-two-samples.md)
3. [Chapter 13: Goodness-of-Fit Tests and Categorical Data Analysis](notes/ch13-goodness-of-fit-categorical-data.md)
4. [Chapter 11: The Analysis of Variance (ANOVA)](notes/ch11-anova.md)
5. [Chapter 12: Regression and Correlation](notes/ch12-simple-linear-regression.md)

All five chapters are also in one PDF: [MATH260_Lecture_Notes.pdf](MATH260_Lecture_Notes_Final.pdf)

## Formats

- `notes/`: the chapters in Markdown. You can read them directly on GitHub.
- `MATH260_Lecture_Notes.pdf`: all chapters in one file.
- `latex/`: the LaTeX source of the PDF. To build it, run `cd latex && latexmk -pdf main.tex`.
- `figures/`: the graphs used in the notes. You can redraw them with `python scripts/make_figs.py` (needs numpy, scipy and matplotlib).