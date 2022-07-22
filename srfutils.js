#!/usr/bin/env k8

/*******************************
 * Command line option parsing *
 *******************************/

var getopt = function(args, ostr) {
	var oli; // option letter list index
	if (typeof(getopt.place) == 'undefined')
		getopt.ind = 0, getopt.arg = null, getopt.place = -1;
	if (getopt.place == -1) { // update scanning pointer
		if (getopt.ind >= args.length || args[getopt.ind].charAt(getopt.place = 0) != '-') {
			getopt.place = -1;
			return null;
		}
		if (getopt.place + 1 < args[getopt.ind].length && args[getopt.ind].charAt(++getopt.place) == '-') { // found "--"
			++getopt.ind;
			getopt.place = -1;
			return null;
		}
	}
	var optopt = args[getopt.ind].charAt(getopt.place++); // character checked for validity
	if (optopt == ':' || (oli = ostr.indexOf(optopt)) < 0) {
		if (optopt == '-') return null; //  if the user didn't specify '-' as an option, assume it means null.
		if (getopt.place < 0) ++getopt.ind;
		return '?';
	}
	if (oli+1 >= ostr.length || ostr.charAt(++oli) != ':') { // don't need argument
		getopt.arg = null;
		if (getopt.place < 0 || getopt.place >= args[getopt.ind].length) ++getopt.ind, getopt.place = -1;
	} else { // need an argument
		if (getopt.place >= 0 && getopt.place < args[getopt.ind].length)
			getopt.arg = args[getopt.ind].substr(getopt.place);
		else if (args.length <= ++getopt.ind) { // no arg
			getopt.place = -1;
			if (ostr.length > 0 && ostr.charAt(0) == ':') return ':';
			return '?';
		} else getopt.arg = args[getopt.ind]; // white space
		getopt.place = -1;
		++getopt.ind;
	}
	return optopt;
}

/************
 * Commands *
 ************/

function srf_cmd_enlong(args) {
	var c, opt = { min_len:200, dbl:false };
	while ((c = getopt(args, "l:d")) != null) {
		if (c == 'l') opt.min_len = parseInt(getopt.arg);
		else if (c == 'd') opt.dbl = true;
	}
	if (args.length == getopt.ind) {
		print("Usage: srfutils.js enlong [options] <circ.fa>");
		print("Options:");
		print("  -l INT    min length [" + opt.min_len + "]");
		print("  -d        double sequences");
		return 1;
	}
	var buf = new Bytes();
	var file = new File(args[getopt.ind]);

	function enlong_seq(seq) {
		if (opt.dbl) return seq + "\n" + seq;
		var t = seq;
		while (t.length < opt.min_len)
			t += seq;
		return t;
	}

	var name = null, comm = null, seq = null;
	while (file.readline(buf) >= 0) {
		var m, l = buf.toString();
		if ((m = /^>(\S+)(.*)/.exec(l)) != null) {
			if (name != null) {
				print(">" + name + comm);
				print(enlong_seq(seq));
			}
			name = m[1], comm = m[2], seq = "";
		} else {
			seq += l;
		}
	}
	if (name != null) {
		print(">" + name + comm);
		print(enlong_seq(seq));
	}
	file.close();
	buf.destroy();
}

function srf_drop_paf(opt, t) {
	for (var i = 1; i <= 3; ++i) t[i] = parseInt(t[i]);
	for (var i = 6; i < t.length; ++i) t[i] = parseInt(t[i]);
	if (t[8] - t[7] < opt.min_len) return true;
	var clip = [t[7], t[6] - t[8]];
	if (t[4] == "+") {
		clip[0] = clip[0] < t[2]? clip[0] : t[2];
		clip[1] = clip[1] < t[1] - t[3]? clip[1] : t[1] - t[3];
	} else {
		clip[1] = clip[1] < t[2]? clip[1] : t[2];
		clip[0] = clip[0] < t[1] - t[3]? clip[0] : t[1] - t[3];
	}
	if (t[8] - t[7] < (t[8] - t[7] + clip[0] + clip[1]) * opt.min_cov) return true;
	if (t[9] < t[10] * opt.min_iden) return true;
	t[11] = (t[10] - t[9] + clip[0] + clip[1]) / (t[10] + clip[0] + clip[1]) * 100;
	return false;
}

function srf_cmd_filter(args) {
	var c, opt = { min_len:190, min_cov:0.9, min_iden:0.0 };
	while ((c = getopt(args, "l:c:d:")) != null) {
		if (c == 'l') opt.min_len = parseInt(getopt.arg);
		else if (c == 'c') opt.min_cov = parseFloat(getopt.arg);
		else if (c == 'd') opt.min_iden = parseFloat(getopt.arg);
	}
	if (args.length == getopt.ind) {
		print("Usage: srfutils.js filter [options] <in.paf>");
		print("Options:");
		print("  -l INT      min alignment length [" + opt.min_len + "]");
		print("  -c FLOAT    min coverage [" + opt.min_cov + "]");
		print("  -d FLOAT    min identity [" + opt.min_iden + "]");
		print("Notes: suggested minimap2 setting:");
		print("  minimap2 -c -N1000000 -f1000 -r100,100 <(srfutils.js enlong circ.fa) reads.fa");
		return 1;
	}
	var buf = new Bytes();
	var file = new File(args[getopt.ind]);
	while (file.readline(buf) >= 0) {
		var line = buf.toString();
		var t = line.split("\t", 12);
		if (srf_drop_paf(opt, t)) continue;
		print(line);
	}
	file.close();
	buf.destroy();
}

function srf_cmd_srtflt(args) {
	var c, opt = { min_len:190, min_cov:0.9, min_iden:0.8, to_count:false, min_cnt:20, thres2nd:0.5 };
	while ((c = getopt(args, "l:c:d:u")) != null) {
		if (c == 'l') opt.min_len = parseInt(getopt.arg);
		else if (c == 'c') opt.min_cov = parseFloat(getopt.arg);
		else if (c == 'd') opt.min_iden = parseFloat(getopt.arg);
		else if (c == 'u') opt.to_count = true;
	}
	if (args.length == getopt.ind) {
		print("Usage: srfutils.js srtflt [options] <in.paf>");
		print("Options:");
		print("  -l INT      min alignment length [" + opt.min_len + "]");
		print("  -c FLOAT    min coverage [" + opt.min_cov + "]");
		print("  -d FLOAT    min identity [" + opt.min_iden + "]");
		print("  -u          count instead of filter");
		print("Notes: suggested minimap2 setting:");
		print("  minimap2 -c -N1000000 -f1000 -r100,100 <(srfutils.js enlong circ.fa) reads.fa | sort -k1,1 -k3,3n");
		return 1;
	}
	var buf = new Bytes();
	var file = args[getopt.ind] == "-"? new File() : new File(args[getopt.ind]);
	var hit = [], a = [];
	while (file.readline(buf) >= 0) {
		var line = buf.toString();
		var t = line.split("\t", 12);
		if (srf_drop_paf(opt, t)) continue;
		t[12] = true;
		// remove out-of-range elements
		while (a.length > 0) {
			if (a[0][0] == t[0] && a[0][3] > t[2])
				break;
			var s = a.shift();
			if (s[12]) {
				if (opt.to_count) hit.push([t[5], t[0], t[11]]);
				else print(s.slice(0, 12).join("\t"));
			}
		}
		// check overlapping records
		for (var i = 0; i < a.length; ++i) {
			var s = a[i];
			if (s[12] == false) continue; // already filtered
			if (s[3] <= t[2]) continue; // no overlap
			if (s[11] <= t[11]) {
				t[12] = false;
				break;
			} else s[12] = false;
		}
		if (t[12]) a.push(t);
	}
	while (a.length > 0) {
		var s = a.shift();
		if (opt.to_count) hit.push([t[5], t[0], t[11]]);
		else print(s.slice(0, 12).join("\t"));
	}
	file.close();
	buf.destroy();

	if (opt.to_count) {
		var a = hit.sort(function(x, y) { return x[0]<y[0]?-1:x[0]>y[0]?1:x[1]<y[1]?-1:x[1]>y[1]?1:0 });
		for (var i = 1, i0 = 0; i <= a.length; ++i) {
			if (i == a.length || a[i][0] != a[i0][0]) {
				var b = [];
				for (var j = i0 + 1, j0 = i0; j <= i; ++j) {
					if (j == i || a[j][1] != a[j0][1]) {
						var avg = 0.0;
						for (var k = j0; k < j; ++k)
							avg += a[k][2];
						avg /= j - j0;
						score = (j - j0) * Math.pow(1 - avg / 100, 4);
						if (j - j0 >= opt.min_cnt)
							b.push([a[j0][1], j - j0, avg, score]);
						j0 = j;
					}
				}
				if (b.length > 0) {
					b = b.sort(function(x, y) { return y[3] - x[3]; });
					var j, thres = b[0][3] * opt.thres2nd;
					for (j = 0; j < b.length; ++j)
						if (b[j][3] < thres) break;
					b.length = j;
					var o = [];
					for (j = 1; j < b.length; ++j)
						o.push(b[j][0] + ":" + b[j][3].toFixed(2));
					print(a[i0][0], b[0][0], b[0][1], b[0][2].toFixed(2), b[0][3].toFixed(2), o.length, o.join("\t"));
					//for (var j = 0; j < b.length; ++j) print(a[i0][0], b[j].join("\t"));
				}
				i0 = i;
			}
			//if (i < a.length) print(a[i].join("\t"));
		}
	}
}

function srf_cmd_paf2bed(args) {
	var c, opt = { min_len:190, min_cov:0.9, min_iden:0.8 };
	while ((c = getopt(args, "l:c:d:p")) != null) {
		if (c == 'l') opt.min_len = parseInt(getopt.arg);
		else if (c == 'c') opt.min_cov = parseFloat(getopt.arg);
		else if (c == 'd') opt.min_iden = parseFloat(getopt.arg);
	}
	if (args.length == getopt.ind) {
		print("Usage: srfutils.js paf2bed [options] <in.paf>");
		print("Options:");
		print("  -l INT      min alignment length [" + opt.min_len + "]");
		print("  -c FLOAT    min coverage [" + opt.min_cov + "]");
		print("  -d FLOAT    min identity [" + opt.min_iden + "]");
		print("Notes: suggested minimap2 setting:");
		print("  minimap2 -c -N1000000 -f1000 -r100,100 <(srfutils.js enlong circ.fa) reads.fa");
		return 1;
	}
	var buf = new Bytes();
	var file = args[getopt.ind] == "-"? new File() : new File(args[getopt.ind]);
	var b = [], sum = {};

	function process(b) {
		b = b.sort(function(x,y) { return x[2]-y[2] });
		var j0 = 0;
		for (var j = 0; j < b.length; ++j) {
			var t = b[j];
			while (j0 < b.length && b[j0][3] <= t[2])
				++j0;
			var st = t[2], en = t[3];
			for (var i = j0; i < j; ++i) {
				var s = b[i];
				if (s[13] <= t[2]) continue;
				if (s[11] < t[11])
					st = st > s[13]? st : s[13];
			}
			for (var i = j0; i < j; ++i) {
				var s = b[i];
				if (s[13] <= st) continue;
				if (s[11] >= t[11])
					s[13] = st;
			}
			t[12] = st, t[13] = en;
		}
		var k = 0;
		for (var j = 0; j < b.length; ++j) {
			var t = b[j];
			if (t[13] - t[12] <= 0) continue;
			b[k++] = t;
		}
		b.length = k;
		if (k == 0) return;
		en = b[0][13];
		var sc = (b[0][13] - b[0][12]) * b[0][11];
		for (var j = 1, j0 = 0; j <= b.length; ++j) {
			if (j == b.length || b[j][0] != b[j0][0] || b[j][5] != b[j0][5] || b[j][12] != en) {
				print(b[j0][0], b[j0][12], en, b[j0][5], sc / (en - b[j0][12]));
				if (j < b.length)
					j0 = j, en = b[j][13], sc = (b[j][13] - b[j][12]) * b[j][11];
			} else {
				en = b[j][13], sc += (b[j][13] - b[j][12]) * b[j][11];
			}
		}
		b.length = 0;
	}

	while (file.readline(buf) >= 0) {
		var line = buf.toString();
		var t = line.split("\t", 12);
		if (srf_drop_paf(opt, t)) continue;
		t[12] = t[2], t[13] = t[3];
		if (b.length > 0 && b[0][0] != t[0])
			process(b);
		b.push(t);
	}
	process(b);
	if (!opt.print)
		for (var x in sum)
			print(x, sum[x]);
	file.close();
	buf.destroy();
}

String.prototype.lpad = function(padString, length) {
	var str = this;
	while (str.length < length)
		str = padString + str;
	return str;
}

function srf_cmd_merge(args) {
	var c, opt = {};
	while ((c = getopt(args, "")) != null) {
	}
	if (args.length - getopt.ind < 2) {
		print("Usage: srfutils.js merge <prefix> <in.fa>");
		return;
	}
	var pre = args[getopt.ind];
	var buf = new Bytes();
	var file = new File(args[getopt.ind + 1]);
	var h = {}, n = 0, seq = '', name = null;
	while (file.readline(buf) >= 0) {
		var m, line = buf.toString();
		if ((m = /^>(\S+)/.exec(line)) != null) {
			if (name != null) {
				if (h[seq] == null) {
					++n;
					print('>' + pre + '*' + String(n).lpad('0', 4) + ' ' + name);
					print(seq);
					h[seq] = 1;
				}
			}
			seq = '', name = m[1];
		} else seq += line;
	}
	if (name != null && h[seq] == null) {
		++n;
		print('>' + pre + '*' + String(n).lpad('0', 4) + ' ' + name);
		print(seq);
	}
	file.close();
	buf.destroy();
}

/*************************
 ***** main function *****
 *************************/

function main(args)
{
	if (args.length == 0) {
		print("Usage: srfutils.js <command> [arguments]");
		print("Commands:");
		print("  enlong    enlong circ contigs");
		print("  filter    filter read-to-circ alignment");
		print("  paf2bed   extract non-overlapping regions with satellites");
		print("  srtflt    filter and count for sorted alignment");
		print("  merge     merge identical sequences");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'filter') srf_cmd_filter(args);
	else if (cmd == 'srtflt') srf_cmd_srtflt(args);
	else if (cmd == 'enlong') srf_cmd_enlong(args);
	else if (cmd == 'merge') srf_cmd_merge(args);
	else if (cmd == 'paf2bed') srf_cmd_paf2bed(args);
	else throw Error("unrecognized command: " + cmd);
}

main(arguments);
