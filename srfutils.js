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
	var iden = t[9] / (t[10] + clip[0] + clip[1]);
	if (iden < opt.min_iden) return true;
	t[11] = (1 - iden) * 100;
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

function srf_cmd_bed2cnt(args) {
	var c, opt = { min_len:10000, thres2nd:0.5 };
	while ((c = getopt(args, "")) != null) {
	}
	if (args.length == getopt.ind) {
		print("Usage: srfutils.js bed2cnt <in.bed>");
		return 1;
	}
	var buf = new Bytes();
	var file = args[getopt.ind] == "-"? new File() : new File(args[getopt.ind]);
	var hit = [];
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		hit.push([t[3], t[0], parseInt(t[2]) - parseInt(t[1]), parseFloat(t[4])]);
	}
	file.close();
	buf.destroy();

	var a = hit.sort(function(x, y) { return x[0]<y[0]?-1:x[0]>y[0]?1:x[1]<y[1]?-1:x[1]>y[1]?1:0 });
	for (var i = 1, i0 = 0; i <= a.length; ++i) {
		if (i == a.length || a[i][0] != a[i0][0]) {
			var b = [];
			for (var j = i0 + 1, j0 = i0; j <= i; ++j) {
				if (j == i || a[j][1] != a[j0][1]) {
					var len = 0, score = 0;
					for (var k = j0; k < j; ++k)
						len += a[k][2], score += a[k][2] * a[k][3];
					var avg = score / len;
					score = len * Math.pow(1 - avg / 100, 4);
					if (len >= opt.min_len)
						b.push([a[j0][1], len, avg, score]);
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
			}
			i0 = i;
		}
	}
}

function srf_cmd_bed2abun(args) {
	var c, opt = { accumu:false, binplot:false, gsize:0, label:"" };
	while ((c = getopt(args, "abg:l:")) != null) {
		if (c == 'a') opt.accumu = true;
		else if (c == 'b') opt.binplot = true;
		else if (c == 'g') opt.gsize = parseInt(getopt.arg);
		else if (c == 'l') opt.label = getopt.arg;
	}
	if (getopt.ind == args.length) {
		print("Usage: srfutils.js bed2abun [options] <in.bed>");
		print("Options:");
		print("  -g INT       total length [0]");
		print("  -a           generate data for accumulative plot");
		print("  -b           for stacked bar plot (<100,<1000,<2000,<5000,<10000,>=10000)");
		print("  -l STR       label with -4");
		return;
	}

	var buf = new Bytes();
	var file = args[getopt.ind] == "-"? new File() : new File(args[getopt.ind]);
	var abun = {};
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (abun[t[3]] == null) abun[t[3]] = [0, 0];
		var len = parseInt(t[2]) - parseInt(t[1]);
		abun[t[3]][0] += len;
		abun[t[3]][1] += len * parseFloat(t[4]);
	}
	file.close();
	buf.destroy();
	var a = [];
	for (var x in abun)
		a.push([x, abun[x][0], abun[x][1] / abun[x][0]]);
	var tot = 0;
	for (var i = 0; i < a.length; ++i)
		tot += a[i][1];
	if (opt.gsize > tot) tot = opt.gsize;
	if (opt.accumu || opt.binplot) {
		var b4 = [0, 0, 0, 0, 0, 0];
		for (var i = 0; i < a.length; ++i) {
			var m;
			if ((m = /-(\d+)$/.exec(a[i][0])) == null)
				throw Error("failed to match format");
			var len = parseInt(m[1]);
			a[i].push(len);
			if (len < 100) b4[0] += a[i][1];
			else if (len < 1000) b4[1] += a[i][1];
			else if (len < 2000) b4[2] += a[i][1];
			else if (len < 5000) b4[3] += a[i][1];
			else if (len < 10000) b4[4] += a[i][1];
			else b4[5] += a[i][1];
		}
		if (opt.accumu) {
			a.sort(function(x,y) { return x[3] - y[3] });
			for (var i = 0, c = 0; i < a.length; ++i) {
				c += a[i][1] / tot;
				print(c, a[i][3], a[i][1]);
			}
		} else {
			var label = opt.label != ""? opt.label : args[getopt.ind];
			for (var i = 0; i < b4.length; ++i)
				b4[i] /= tot;
			print(label, b4.join("\t"));
		}
	} else {
		a = a.sort(function(x,y) { return y[1] - x[1] });
		for (var i = 0; i < a.length; ++i)
			print(a[i].join("\t"), a[i][1] / tot);
	}
}

function srf_cmd_paf2bed(args) {
	var c, opt = { min_len:190, min_cov:0.9, min_iden:0.9 };
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
		print("  bed2cnt   per-chromosome count");
		print("  bed2abun  per-contig abundance");
		print("  merge     merge identical sequences");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'filter') srf_cmd_filter(args);
	else if (cmd == 'enlong') srf_cmd_enlong(args);
	else if (cmd == 'merge') srf_cmd_merge(args);
	else if (cmd == 'paf2bed') srf_cmd_paf2bed(args);
	else if (cmd == 'bed2cnt') srf_cmd_bed2cnt(args);
	else if (cmd == 'bed2abun') srf_cmd_bed2abun(args);
	else throw Error("unrecognized command: " + cmd);
}

main(arguments);
