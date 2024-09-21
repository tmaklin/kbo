#[test]
fn map_nissle_against_clbs() {
    let (sbwt, lcs) = sablast::index::build_sbwt(&"tests/data/clbS.fna.gz".to_string(), &None);

    let expected = vec![(455, 967, 512, 1, '+'),(996, 1001, 4, 2, '+'),(998, 1000, 3, 0, '-')];
    let got = sablast::map(&"tests/data/NZ_CP058217.1_clbS.fna.gz".to_string(), &sbwt::SbwtIndexVariant::SubsetMatrix(sbwt), &lcs.unwrap());

    assert_eq!(got, expected);
}
