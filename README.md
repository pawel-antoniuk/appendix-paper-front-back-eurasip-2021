# (Appendix) Automatic discrimination between front and back ensemble locations in HRTF-convolved binaural recordings of music 
This repository consists of scripts and data useful to replicate the experiments described in the paper "Automatic discrimination between front and back ensemble locations in HRTF-convolved binaural recordings of music" (EURASIP Journal on Audio, Speech, and Music Processing 2021).

## Structure
The repository is organized as follows:
- [classification-traditional](classification-traditional) contains scripts and data used in the development of traditional classification algorithms as described in the section "5.1. Traditional classification algorithms" of the previously mentioned article. 
- [classification-deep](classification-deep) contains sources and data used in the development of the deep learning algorithm as described in the section "5.2. Deep learning algorithm".
- [convolver](convolver) contains the source code responsible for generating binaural audio files used in the experiments.

## Citation
If you find this repository useful for your research, please consider citing our paper using the "Cite this repository" option on the repository's GitHub page or using the following BibTeX entry:

```bibtex
@article{Zielinski_Antoniuk_Hyunkook_Dale_2021, 
    title={Automatic discrimination between front and back ensemble locations in HRTF-convolved binaural recordings of music}, 
    journal={EURASIP Journal on Audio, Speech, and Music Processing}, 
    author={Zieliński, Sławomir K. and Antoniuk, Paweł and Hyunkook, Lee and Dale, Johnson}, 
    year={2021}
}
```

## Dependencies
Software dependencies:
- [Python 3.6+](https://docs.python.org/3.6/) - a development environment used to implement the traditional algorithm
- [Auditory front-end](http://docs.twoears.eu/en/1.5/afe/) - a software used to extract features that were then given to the traditional algoritms
- [scikit-learn](https://scikit-learn.org/stable/) - a machine learning library used to implement traditional algorithms
- [MATLAB](https://www.mathworks.com/products/matlab.html) - a development environment used to implement the deep learning algorithm
- [VOICEBOX](http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html) - a toolbox used to implement the binaural convolver

## Authors
Sławomir K. Zieliński <sup>1</sup>, Paweł Antoniuk <sup>1</sup>, Hyunkook Lee <sup>2</sup>, and Dale Johnson <sup>2</sup>

<sup>1</sup> Faculty of Computer Science, Białystok University of Technology, 15-351 Białystok, Poland; s.zielinski@pb.edu.pl (S.K.Z.); p.antoniuk6@student.pb.edu.pl (P.A.)

<sup>2</sup> Applied Psychoacoustics Laboratory (APL), University of Huddersfield, Huddersfield HD1 3DH, UK; H.Lee@hud.ac.uk (H.L.); D.S.Johnson2@hud.ac.uk (D.J.)

## License
The content of this repository is licensed under the terms of the GNU General Public License v3.0 license. Please see the [LICENSE](LICENSE) file for more details.
