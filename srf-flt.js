#!/usr/bin/env k8

function main(args) {
	var opt = { min_len:200, min_cov:0.9, min_iden:0.9 };
	var buf = new Bytes();
	var file = args.length > 0? new File(args[0]) : new File();
	while (file.readline(buf) >= 0) {
		var line = buf.toString();
		var t = line.split("\t", 12);
		for (var i = 1; i <= 3; ++i) t[i] = parseInt(t[i]);
		for (var i = 6; i < t.length; ++i) t[i] = parseInt(t[i]);
		if (t[8] - t[7] < opt.min_len) continue;
		var clip = [t[7], t[6] - t[8]];
		if (t[4] == "+") {
			clip[0] = clip[0] < t[2]? clip[0] : t[2];
			clip[1] = clip[1] < t[1] - t[3]? clip[1] : t[1] - t[3];
		} else {
			clip[1] = clip[1] < t[2]? clip[1] : t[2];
			clip[0] = clip[0] < t[1] - t[3]? clip[0] : t[1] - t[3];
		}
		if (t[8] - t[7] < (t[8] - t[7] + clip[0] + clip[1]) * opt.min_cov) continue;
		if (t[9] < t[10] * opt.min_iden) continue;
		print(line);
	}
	file.close();
	buf.destroy();
}

main(arguments);
