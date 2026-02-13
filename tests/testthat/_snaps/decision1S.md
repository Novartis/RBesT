# decision1S works for lower sided

    Code
      print(dec)
    Output
      1 sample decision function
      Conditions for acceptance:
      P(theta <= 0) > 0.05
      P(theta <= 0.8) > 0.5
      P(theta <= 1.2) > 0.05

---

    list(lower = c(2.30258509299405, 0.00636272327712173, 2.31211393398714
    ))

# decision1S works for upper sided

    Code
      print(dec)
    Output
      1 sample decision function
      Conditions for acceptance:
      P(theta > 0) > 0.05
      P(theta > 0.8) > 0.5
      P(theta > 1.2) > 0.05

---

    list(upper = c(2.30258509299405, -0.00640346690337834, 2.29296457895205
    ))

# decision1S works for two sided

    Code
      print(decMixed)
    Output
      1 sample decision function (two-sided)
      Conditions for acceptance:
      Lower tail conditions:
      P(theta <= 0) > 0.05
      P(theta <= 1.2) > 0.05
      Upper tail conditions:
      P(theta > 0.8) > 0.5

---

    list(lower = c(2.30258509299405, 2.31211393398714), upper = -0.00640346690337834)

