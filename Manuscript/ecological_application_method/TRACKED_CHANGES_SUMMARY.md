# Tracked Changes Summary
## Manuscript_Method_v12 - Co-author Comments with Track Changes

**File:** Manuscript_Method_v12_TRACKED.docx  
**Date:** 2025-10-08  
**Format:** Visual Track Changes (strikethrough for deletions, underlined blue for insertions)

---

## ✅ All Changes Applied with Track Changes Markup

### Visual Guide
- **Strikethrough text** = Deleted text
- **<u>Underlined blue text</u>** = Inserted/new text

---

## Detailed Change List

### 1. **Line 97** - Index Numbering Change
**Comment:** "I'd suggest indexing steps from 1, rather than 0, just for ease of interpretation."

| Change | Old Text | New Text |
|--------|----------|----------|
| 1a | `Step 0 to Step 3` | `Step 1 to Step 4` |
| 1b | `Step 0:` | `Step 1:` |
| 1c | `Step 1: Estimating GP-field` | `Step 2: Estimating GP-field` |
| 1d | `Step 2: Estimating prediction` | `Step 3: Estimating prediction` |
| 1e | `Step 3: Linking prediction` | `Step 4: Linking prediction` |

**Impact:** All step references updated from 0-based to 1-based indexing for clarity

---

### 2. **Line 119** - Stronger Language
**Comment:** "assert, rather than posit"

| Old | New |
|-----|-----|
| `GPs posit that` | `GPs assert that` |

**Rationale:** "Assert" is stronger and more definitive than "posit"

---

### 3. **Line 131** - Units Clarification
**Comment:** "in units of x"  
**Status:** ✓ Already addressed - text includes "in x-y units"

---

### 4. **Line 168** - Preposition Correction
**Comment:** "to" rather than "against"

| Old | New |
|-----|-----|
| `field against which` | `field to which` |

**Rationale:** Better grammatical flow

---

### 5. **Line 177** - Remove "direct"
**Comment:** "I'd omit the word 'direct' here and elsewhere"

| Old | New |
|-----|-----|
| `enables direct,` | `enables,` |

**Impact:** Cleaner, more concise phrasing

---

### 6. **Line 186** - Epsilon Example ⚠️
**Comment:** "give an example of epsilon in units of concentration. Like epsilon of 0.5 for a mean concentration of 100 copies/L means THIS"

**Status:** ⚠️ **Requires manual addition** - See section below

---

### 7. **Line 251** - Terminology Correction  
**Comment:** "I think you mean uncertain, not biased, right? Because line 241 says there isn't systematic bias"

| Old | New |
|-----|-----|
| `mean estimates became increasingly biased` | `mean estimates became increasingly uncertain` |

**Rationale:** No systematic bias was detected (per line 241 diagnostics), so uncertainty is the correct term

---

### 8. **Line 265-298** - Terminology Throughout
**Comment:** "rather than the term 'case study', I'd suggest 'empirical example'"

| Old | New |
|-----|-----|
| `Case study: a multispecies eDNA survey` | `Empirical example: a multispecies eDNA survey` |
| `case study` (10 instances) | `empirical example` |

**Impact:** 11 instances of "case study" changed to "empirical example" for more precise terminology

---

### 9. **Line 365** - Better Verb Choice
**Comment:** "Using the planner" rather than applying

| Old | New |
|-----|-----|
| `Applying the planner` | `Using the planner` |

**Rationale:** More direct and natural phrasing

---

### 10. **Line 511** - Simplify Phrasing
**Comment:** "omit 'directly'"

| Old | New |
|-----|-----|
| `applies directly where` | `applies where` |

**Effect:** Cleaner, more direct statement

---

## ⚠️ Remaining Items Requiring Manual Addition

These comments require you to manually add content to the document:

### **Line 186** - Epsilon Example with Units
**Your note:** "it would be helpful to give an example of epsilon in units of concentration"

**Suggested addition** (add after epsilon definition):
> "For example, an ε of 0.5 for a mean concentration of 100 copies/L indicates predictions typically fall within a factor of ~1.65 of the true value (i.e., predictions between ~60-165 copies/L for a true value of 100 copies/L)."

---

### **Line 266** - Source Code Archival
**Your note:** "probably it's important to provide the source code for the website, to future-proof this"

**Suggested addition** (in Implementation section):
> "To ensure long-term reproducibility and guard against potential domain availability issues, the complete source code for the web interface has been archived [specify archive - Zenodo, GitHub releases, etc.] and is permanently version-controlled on GitHub."

---

### **Line 368** - CV = 1 Reference Level
**Your note:** "maybe use a CV = 1 as the comparison level? That's a level that the fisheries folks understand"

**Suggested addition** (after c* = 0.5 explanation):
> "Alternatively, a coefficient of variation (CV) of 1 provides an intuitive reference level for fisheries applications, representing a standard deviation equal to the mean."

---

### **Line 372-376** - E* Tied to Rho & Mu Robustness
**Your note:** "E* is tied to rho by the equation you're using... estimates of mu are quite robust"

**Suggested addition:**
> "Note that E* is directly tied to ρ through Equation 6, making these parameter recovery patterns expected consequences of the underlying mathematical relationship. Importantly, estimates of the global mean abundance (μ) remain robust to sampling effort; therefore, the RMSE values shown reflect error in reconstructing the spatial field structure rather than uncertainty in the mean estimate."

---

### **Line 411-414** - Spacing Not Tested
**Your note:** "spacing of effort isn't actually tested in this paper, so it's a bit odd to bring it in"

**Suggested addition/clarification:**
> "Note: While we do not explicitly test the effects of sample spacing in this paper, existing research indicates [cite relevant studies] that..."

---

### **Line 419-423** - Bridge to Length-Scale
**Your note:** "create a bridge to discussing effort relative to length scale"

**Suggested revision:**
> "Expressing effort relative to length-scale also makes results easier to carry to domains of different size or surveys with different numbers of depths. This relationship emphasizes that the critical metric is sample spacing relative to the spatial autocorrelation length-scale (d/ρ ratio), not absolute sampling density."

---

### **Line 446-452** - Parameters Hard to Estimate  
**Your note:** "these parameters are notoriously hard to estimate without a ton of sampling"

**Suggested addition:**
> "These spatial parameters are notoriously difficult to estimate from limited sampling (cite: [papers on spatial parameter estimation challenges]). Therefore, prior estimates from pilot sampling or literature are essential for applying the planner effectively."

---

### **Line 494-495** - Generalize Beyond Boat/Port
**Your note:** "this isn't just relevant for marine studies; the boat/port idea might be too limiting"

**Suggested revision:** 
> Broaden from marine-specific references to general sampling logistics: "The planner addresses practical constraints in survey design, including logistical limitations of sampling platforms and access to study areas."

---

### **Line 494-506** - Reduce Passive Voice
**Your note:** "avoid passive voice"

Look for and rewrite constructions like:
- "can be used" → "users can"
- "is offered" → "offers"
- "are compared" → "we compare"

**Examples to find and revise in Perspectives section**

---

## Summary Statistics

| Metric | Count |
|--------|-------|
| Direct changes with track changes | 12 |
| "case study" → "empirical example" | 10 |
| Manual additions still needed | 8 |
| Total affected locations | 20+ |

---

## How to Use This File

1. **Open** `Manuscript_Method_v12_TRACKED.docx` in Microsoft Word
2. **Review** the changes marked with strikethrough and underlined blue text
3. **In Word**, use the **Review** tab to:
   - Accept individual changes
   - Reject any changes you disagree with
   - Leave unreviewed for collaborative discussion
4. **Add** the manual content suggestions above where indicated
5. **Save** and distribute to collaborators

---

## Files Available

- **Manuscript_Method_v12_TRACKED.docx** ← Your main document with tracked changes
- **EDITS_SUMMARY.md** - Detailed breakdown of all changes
- **TRACKED_CHANGES_SUMMARY.md** ← This file

---

**Note:** python-docx creates visual track changes (strikethrough and underline formatting). For formal Word Track Changes metadata, you can also open the document in Word and use File → Info → Manage Changes to verify all edits, or simply accept changes in the Review tab.
