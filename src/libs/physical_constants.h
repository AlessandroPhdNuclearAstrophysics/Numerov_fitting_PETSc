  static const double hbar = 197.327;
  static const double Mp = 938.272, Mn = 939.565, Mmu = Mp*Mn/(Mn+Mp);
  static const double htm = hbar*hbar/(2*Mmu);

  double k2_from_E(const double E) {
    double k2;
    k2 = E/htm;
    return k2;
  }