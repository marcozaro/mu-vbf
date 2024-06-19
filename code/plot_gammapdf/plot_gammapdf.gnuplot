set output "gampdf.pdf"
plot \
"gampdf_mmu.txt" u 1:2 t "{\Symbol G}_1, Q={\Symbol m}",\
"gampdf_mmu.txt" u 1:3 t "IWW Q={\Symbol m}",\
"gampdf_mmu.txt" u 1:4 t "NNL-MSb Q={\Symbol m}"
