fn main() {
    cc::Build::new()
        .file("../power_sums.c")
        .compile("power_sums");
}
