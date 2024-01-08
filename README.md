# B-SIDH-attack

・B-SIDH attackをしたい場合.

load_this.m をloadすると動きます. (こちらがアーベル多様体(遅い))

or

load_this_new.mをloadすると動きます. (こちらがKummer(速い)).


・p,N_A,N_Bの具体例について.

example_parameter.mの中にp,N_A,N_Bの例が入っています.


・素数次数同種写像の計算時間を測りたい場合.

load_this.mの中の上から6個(//----より上まで)をloadして下さい. その後, test_isogeny_time.m内の関数compute_isogeny(p,l)を用いてください. test_isogeny_time.mの下あたりに実装例を書いています.
