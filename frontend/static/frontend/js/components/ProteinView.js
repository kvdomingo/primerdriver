import React, { Component } from "react";
import { withRouter } from "react-router-dom";
import LoadingScreen from "./shared/LoadingScreen";
import { Form, ProteinSequenceInput, MutationSelector, MutationType, AdvancedSettings } from "./form";

class ProteinView extends Component {
  state = {
    cursorPosition: 1,
    sequenceLength: 0,
    isValid: false,
    loading: false,
    mode: "PRO",
    sequence: "",
    mutation_type: "",
    target: "",
    position: 0,
    replacement: "",
    Tm_range_min: 75,
    Tm_range_max: 85,
    gc_range_min: 40,
    gc_range_max: 60,
    length_min: 25,
    length_max: 45,
    flank5_range_min: 11,
    flank5_range_max: 21,
    flank3_range_min: 11,
    flank3_range_max: 21,
    forward_overlap5: 9,
    forward_overlap3: 9,
    terminate_gc: true,
    center_mutation: true,
    primer_mode: "complementary",
    expression_system: "Homo sapiens",
    expressionList: [],
  };

  formDefaults = { ...this.state };

  componentDidMount() {
    fetch("/api/expressionsys")
      .then(res => res.json())
      .then(res => {
        let expressionList = res.data;
        this.setState({ expressionList });
      });
  }

  handleChange = e => {
    let { name, value } = e.target;
    if (name === "sequence" || name === "target" || name === "replacement") {
      let sequence = value.toUpperCase().split("");
      let aminoAcids = [
        "A",
        "R",
        "N",
        "D",
        "B",
        "C",
        "E",
        "Q",
        "Z",
        "G",
        "H",
        "I",
        "L",
        "K",
        "M",
        "F",
        "P",
        "S",
        "T",
        "W",
        "Y",
        "V",
      ];
      let filteredSequence = [];
      sequence.forEach(char => {
        if (aminoAcids.includes(char)) filteredSequence.push(char);
      });
      if (e.target.name === "sequence") {
        this.setState({
          sequence: filteredSequence.join(""),
          cursorPosition: e.target.selectionStart + 1,
          sequenceLength: value.length,
        });
      } else {
        this.setState({ [name]: filteredSequence.join("") });
      }
    } else {
      if (name === "mutation_type") {
        this.setState({
          target: "",
          replacement: "",
        });
      }
      this.setState({ [name]: value });
    }
  };

  handleChangeInt = e => {
    let { name, value } = e.target;
    this.setState({ [name]: parseInt(value) });
  };

  handleChangeFloat = e => {
    let { name, value } = e.target;

    if (name === "Tm_range_min" && value >= this.state.Tm_range_max) return;
    if (name === "Tm_range_max" && value <= this.state.Tm_range_min) return;
    if (name === "gc_range_min" && value >= this.state.gc_range_max) return;
    if (name === "gc_range_max" && value <= this.state.gc_range_min) return;

    this.setState({ [name]: parseFloat(value) });
  };

  handleChangeChecked = e => {
    let { name, checked } = e.target;
    this.setState({ [name]: checked });
  };

  handleReset = () => {
    this.setState({ ...this.formDefaults });
  };

  handleSubmit = e => {
    this.setState({ loading: true });
    const xhr = new XMLHttpRequest();
    xhr.open("POST", "/api");
    xhr.onload = () => {
      let res;
      if (xhr.status === 200) {
        res = JSON.parse(xhr.responseText);
      } else {
        res = "Request failed. Please try again later.";
      }
      this.props.responseCatcher(res, this.state.mode);
      let { history } = this.props;
      history.push("/results");
    };
    const data = new FormData();
    Object.keys(this.formData).map(key => data.append(key, this.formData[key]));
    data.append("settings", JSON.stringify(this.formSettings));
    xhr.send(data);
    e.preventDefault();
  };

  validateForm = e => {
    let validSequence = this.state.sequenceLength > 0;
    let validMutation = this.state.mutation_type !== "";
    let validMutationCode =
      this.state.mutation_type === "sub"
        ? this.state.target.length > 0 && this.state.replacement.length > 0
        : this.state.target.length > 0 || this.state.replacement.length > 0;
    validSequence && validMutation && validMutationCode
      ? this.setState({ isValid: true })
      : this.setState({ isValid: false });
    let formData = ["mode", "sequence", "target", "position", "replacement", "mutation_type"];
    this.formData = {};
    formData.map(item => (this.formData[item] = this.state[item]));
    let advSettings = [
      "Tm_range_min",
      "Tm_range_max",
      "gc_range_min",
      "gc_range_max",
      "length_min",
      "length_max",
      "flank5_range_min",
      "flank5_range_max",
      "flank3_range_min",
      "flank3_range_max",
      "forward_overlap5",
      "forward_overlap3",
      "terminate_gc",
      "center_mutation",
      "primer_mode",
      "expression_system",
    ];
    this.formSettings = {};
    advSettings.map(item => (this.formSettings[item] = this.state[item]));
  };

  render() {
    if (this.state.loading) return <LoadingScreen />;
    else
      return (
        <Form
          title="Protein-based Primer Design"
          handleValidate={this.validateForm}
          handleSubmit={this.handleSubmit}
          handleReset={this.handleReset}
          {...this.state}
        >
          <ProteinSequenceInput handleChange={this.handleChange} {...this.state} />
          <MutationType handleChange={this.handleChange} {...this.state} />
          <MutationSelector handleChange={this.handleChange} {...this.state} />
          <AdvancedSettings
            handleChange={this.handleChange}
            handleChangeInt={this.handleChangeInt}
            handleChangeFloat={this.handleChangeFloat}
            handleChangeChecked={this.handleChangeChecked}
            {...this.state}
          />
        </Form>
      );
  }
}

export default withRouter(ProteinView);
