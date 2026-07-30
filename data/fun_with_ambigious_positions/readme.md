## fun with ambigious positions

These files are (very) loosely based on an L7 sample from SRA plus some manually masked sites. They demonstrate how the way UShER + matUtils handle ambigious positions sometimes appears inconsistent. Names are based on the phonetic alphabet.

### India and Juliett: Masked positions meet expectations
India and Juliett are truncated versions of Charlie and Delta. India consists of single SNP, and the same seven masked positions. Juliett shares India's SNP but not the masked positions.

```
>INDIA
C	1960	1
-	1977	1
-	2532	1
-	4013	1
-	7362	1
-	7585	1
-	8876	1
-	9143	1
>JULIETT
C	1960	1
```

What you would expect is that India and Juliett to be 0 apart (or perhaps 7 apart), and after placing them on the 70-sample tree and running matOptimize, they will indeed be zero apart.

### Charlie and Delta: Reversions to reference
Charlie and Delta are identical except for that Charlie has the same seven masked positions India has, while Delta does not. They are much longer than India and Juliett, but you might think that doesn't matter since all that extra stuff is identical between Charlie and Delta. Alas, it does matter.

When placed on the 70-sample tree, Charlie and Delta will consistently be placed 6 apart, before and after matOptimize.

The reason for this is that Charlie and Delta's other positions, which are not present in India and Juliett, place them at a place on the tree that has mutations at positions such as A1977G. However, Charlie asserts a mask at 1977 and Delta implicitly calls reference at 1977, which means Delta has a reversion to reference (G1977A) that Charlie lacks.

To quote Angie Heinrichs: "A sample with masked positions is placed based on its non-masked substitutions, with masked positions inferred."

## Conclusion
A masked position is assumed to be whatever would be expected from previous nodes.
	* In the case of India and Juliett, there is no previous node with a mutation at 1977, so when UShER sees India calling mask at 1977 and Juliett calling reference at 1977, it assumes both are reference.
	* In the case of Charlie and Delta, there is a previous node with A1977G, so when UShER sees a masked position in Charlie, it assumes Charlie is also A1977G. On the other hand, Delta implicitly calls reference.