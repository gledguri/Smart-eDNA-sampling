# Manuscript Edits Summary
## Manuscript_Method_v12 - Co-author Comments Incorporated

**Date:** 2025-10-08  
**Edited by:** Claude with python-docx  
**New file:** Manuscript_Method_v12_REVISED.docx

---

## ✅ Automated Edits Completed

### Successfully Replaced:

1. **Line 97** ✓ 
   - Changed: "Step 0 to Step 3" → "Step 1 to Step 4" (1 instance)
   - Changed: "Step 0:" → "Step 1:" (1 instance)
   - *Effect: Updated indexing from 0-based to 1-based for clarity*

2. **Line 119** ✓
   - Changed: "GPs posit" → "GPs assert"
   - *Effect: Stronger language per co-author preference*

3. **Line 168** ✓
   - Changed: "field against which" → "field to which"
   - *Effect: Improved preposition usage*

4. **Line 177** ✓
   - Removed: "direct" from "enables direct"
   - Changed: "enables direct," → "enables,"
   - *Effect: Cleaner phrasing*

5. **Line 251** ✓
   - Changed: "biased" → "uncertain" (regarding mean estimates)
   - *Effect: More accurate terminology (no systematic bias per line 241 checks)*

6. **Line 265-298** ✓
   - Changed: "Case study:" → "Empirical example:" (2 instances)
   - Changed: "case study" → "empirical example" (8 additional instances throughout)
   - *Effect: More precise terminology for study application*

7. **Line 365** ✓
   - Changed: "Applying the planner" → "Using the planner"
   - *Effect: Improved verb choice*

8. **Line 511** ✓
   - Removed: "directly" from "applies directly where"
   - Changed: "applies directly where" → "applies where"
   - *Effect: Simplified phrasing*

---

## ⚠️ Remaining Items Requiring Manual Review/Addition

### 1. **Line 131** - Units Clarification
**Comment:** "in units of x"  
**Current text:** "ρ is the length-scale parameter (in x-y units; hereafter"  
**Status:** Already adequately addresses the comment with "in x-y units"  
**Action:** CONFIRM - text appears sufficient as-is

### 2. **Line 186** - Epsilon Example ⚠️ NEEDS ADDITION
**Comment:** "I think it would be helpful to give an example of epsilon in units of concentration. Like epsilon of 0.5 for a mean concentration of 100 copies/L means THIS, etc."  
**Current text:** Section defines ε in prediction error equation  
**Action NEEDED:** 
- Add concrete example after the epsilon definition
- Suggested addition: "For example, an ε of 0.5 for a mean concentration of 100 copies/L indicates that predictions typically fall within a factor of ~1.65 of the reference value (100 copies/L will be predicted between 60-165 copies/L)."

### 3. **Line 251-252** - Spatial Structure Emphasis
**Comment:** "'I think it would be helpful to give an example of epsilon in units of concentration. Like epsilon of 0.5 for a mean concentration of 100 copies/L means THIS, etc.'"  
**Current changes:** "biased" → "uncertain" ✓  
**Additional note:** Comment mentions "LESS spatially structured, not more — as rho increases, things get more smooth, and so less spatially structured"  
**Action NEEDED:** 
- Review context around parameter recovery section
- Ensure clarity that higher ρ = smoother/less spatially structured fields

### 4. **Line 266** - Source Code Archival ⚠️ NEEDS ADDITION  
**Comment:** "probably it's important to provide the source code for the website, to future-proof this. Like, what happens when the hosting domain doesn't host it anymore, or something breaks in the javascript, etc?"  
**Current text:** Implementation section describes the web-based tool  
**Action NEEDED:**
- Add paragraph about source code archival strategy
- Should mention: GitHub repository, Zenodo archival (already mentioned in open research statement?)
- Suggested text: "To ensure long-term reproducibility, the complete source code for the eDNA Survey Planner web interface has been archived in [archive service] and remains version-controlled on GitHub. This archival strategy mitigates risks from domain name changes or future hosting unavailability."

### 5. **Line 368** - CV = 1 Comparison Level ⚠️ NEEDS CONSIDERATION
**Comment:** "maybe use a CV = 1 as the comparison level? That's a level that the fisheries folks understand pretty intuitively."  
**Current text:** Shows c* = 0.5 as "arbitrary value chosen for illustration"  
**Action NEEDED:**
- Consider adding note that CV = 1 represents a reference level fisheries managers understand (e.g., coefficient of variation = 1 means std dev = mean)
- Could be added as parenthetical or separate sentence

### 6. **Line 372-376** - E* Tied to Rho Relationship ⚠️ NEEDS EMPHASIS
**Comment:** "E* is tied to rho by the equation you're using, so these are all quite expected — it might be worth stating that explicitly. Moreover, it's probably worth emphasizing — as you do elsewhere — that estimates of mu are quite robust to sampling effort. So while the RMSE here reflects overall error about the spatial field and distribution of concentration, it could be misleading if people read it as a measure of overall error (i.e., in the mean, which is often an important parameter of interest)."  
**Action NEEDED:**
- Add explicit statement that E* is a function of ρ (from Equation)
- Emphasize that μ estimates are robust to sampling effort
- Clarify that RMSE measures field reconstruction error, not mean estimation error
- Suggested addition: "Note that E* is directly tied to ρ through Equation X, making these results expected consequences of the mathematical relationship. Importantly, estimates of the global mean abundance (μ) remain robust to sampling effort, so the RMSE values here reflect error in reconstructing the spatial field rather than uncertainty in the mean estimate."

### 7. **Line 411-414** - Spacing Not Tested ⚠️ NEEDS CLARIFICATION
**Comment:** "spacing of effort isn't actually tested in this paper, so it's a bit odd to bring it in here."  
**Action NEEDED:**
- If spacing is discussed but not tested, add clarifying note
- Suggested text: "While we do not explicitly test the effects of sample spacing in this paper, research indicates that..."

### 8. **Line 419-423** - Bridge to Length-Scale Discussion ⚠️ NEEDS ENHANCEMENT
**Comment:** "same thing here; create a bridge to discussing effort relative to length scale, so it becomes a real point of discussion, rather than what might seem like a sneaky result."  
**Action NEEDED:**
- Create clearer transition between effort discussion and length-scale concepts
- Emphasize the relationship: effort should be discussed relative to ρ, not in isolation
- Consider adding discussion of d/ρ ratio (sample spacing to length-scale)

### 9. **Line 446-452** - Parameter Estimation Difficulty ⚠️ NEEDS CITATIONS
**Comment:** "these parameters are, I believe, notoriously hard to estimate without a ton of sampling, and you can probably cite literature for that idea."  
**Action NEEDED:**
- Add citations about difficulty of estimating spatial parameters
- Suggested literature topics: spatial autocorrelation estimation challenges, small sample properties of variogram/correlation estimators
- Could add: "These parameters are notoriously difficult to estimate without substantial sampling effort (cite papers on spatial parameter estimation with small samples)."

### 10. **Line 494-495** - Generalize Beyond Boat/Port ⚠️ NEEDS BROADENING
**Comment:** "note this isn't just relevant for marine studies; the boat/port idea might be too limiting."  
**Action NEEDED:**
- Broaden language beyond marine/boat/port context
- Should address sampling logistics in general
- Suggested revision: Consider if point about surveys/planning applies to terrestrial, aquatic, or other contexts
- Perhaps replace "boat/port" with "sampling platform/logistics" or similar

### 11. **Line 494-506** - Passive Voice Reduction ⚠️ NEEDS REVIEW
**Comment:** "avoid passive voice" throughout this section  
**Action NEEDED:**
- Review "Perspectives and applications" section (lines 494-506 and beyond)
- Convert passive constructions to active voice
- Examples to look for:
  - "can be used" → "users can"
  - "is offered" → "offers"  
  - "are intended" → "intends to"
  - Etc.
- This requires careful reading of the full section

---

## Summary Statistics

| Category | Count |
|----------|-------|
| Direct text replacements completed | 10 |
| Section changes completed | 1 |
| Instances of "case study" → "empirical example" | 10 |
| Items requiring manual review | 11 |
| Items requiring manual addition | 6 |

---

## Recommended Next Steps

1. **High Priority** (clarifying edits):
   - [ ] Add epsilon example with units (Line 186)
   - [ ] Add source code archival note (Line 266)
   - [ ] Add E* = f(ρ) and μ robustness clarifications (Line 372-376)

2. **Medium Priority** (enhancing clarity):
   - [ ] Review/emphasize spatial structure discussion (Line 251)
   - [ ] Add clarification on spacing not tested (Line 411)
   - [ ] Create better effort/length-scale bridge (Line 419-423)
   - [ ] Generalize beyond boat/port (Line 494-495)

3. **Lower Priority** (polishing):
   - [ ] Add citations on parameter estimation difficulty (Line 446-452)
   - [ ] Consider CV = 1 note for fisheries context (Line 368)
   - [ ] Review and reduce passive voice (Line 494-506)
   - [ ] Confirm units clarification adequate (Line 131)

---

## Files

- **Original:** Manuscript_Method_v12.docx
- **Revised:** Manuscript_Method_v12_REVISED.docx (all automated edits applied)
- **PDF Reference:** Manuscript_Method_v12.pdf

---

**Note:** The revised document is ready for the remaining manual edits. All straightforward text replacements have been applied automatically.
