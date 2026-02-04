# decision2S works for lower sided

    Code
      print(decLower)
    Output
      2 sample decision function
      Conditions for acceptance:
      P(theta1 - theta2 <= 0) > 0.05
      P(theta1 - theta2 <= 0.8) > 0.5
      Link: identity 

---

    list(lower = c(2.31037842239559, 0.0139840577654454))

# decision2S works for upper sided

    Code
      print(decUpper)
    Output
      2 sample decision function
      Conditions for acceptance:
      P(theta1 - theta2 > 1.2) > 0.05
      P(theta1 - theta2 > 2) > 0.5
      Link: identity 

---

    list(upper = c(2.28522400871107, -0.0237480873584778))

# decision2S works for two sided

    Code
      print(decMixed)
    Output
      2 sample two-sided decision function
      Lower side conditions for acceptance:
      P(theta1 - theta2 <= 0) > 0.05
      P(theta1 - theta2 <= 1.2) > 0.05
      Upper side conditions for acceptance:
      P(theta1 - theta2 > 0.8) > 0.5
      Link: identity 

---

    list(lower = c(2.31037842239559, 2.31964990630052), upper = -0.0141823883816593)

