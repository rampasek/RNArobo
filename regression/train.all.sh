for file in *.data.out; do
	echo "t = read.table(\"$file\"); fit <- lm(V4 ~ V2 + V3, data=t); coefficients(fit);" | R --vanilla --silent | tail -2 | head -1
done
