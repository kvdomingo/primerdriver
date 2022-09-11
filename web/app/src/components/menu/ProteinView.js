import { useEffect, useState } from "react";
import { useNavigate } from "react-router-dom";
import PropTypes from "prop-types";
import api from "../../api";
import { usePrimerDriverContext } from "../../contexts/PrimerDriverContext";
import { AdvancedSettings, Form, MutationSelector, MutationType, ProteinSequenceInput } from "../form";
import LoadingScreen from "../shared/LoadingScreen";

const AMINO_ACIDS = [
  "A",
  "R",
  "N",
  "D",
  "C",
  "Q",
  "E",
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

function ProteinView(props) {
  const [formData, setFormData] = useState({});
  const [cursorPosition, setCursorPosition] = useState(1);
  const [sequenceLength, setSequenceLength] = useState(0);
  const [isValid, setIsValid] = useState(false);
  const [loading, setLoading] = useState(false);
  const [sequence, setSequence] = useState("");
  const [mutationType, setMutationType] = useState("");
  const [target, setTarget] = useState("");
  const [position, setPosition] = useState(0);
  const [replacement, setReplacement] = useState("");
  const [gcRange, setGcRange] = useState([40, 60]);
  const [tmRange, setTmRange] = useState([75, 85]);
  const [flank5Range, setFlank5Range] = useState([11, 21]);
  const [flank3Range, setFlank3Range] = useState([11, 21]);
  const [lengthRange, setLengthRange] = useState([25, 45]);
  const [forwardOverlap5, setForwardOverlap5] = useState(9);
  const [forwardOverlap3, setForwardOverlap3] = useState(9);
  const [terminateGc, setTerminateGc] = useState(true);
  const [centerMutation, setCenterMutation] = useState(true);
  const [primerMode, setPrimerMode] = useState("complementary");
  const [expressionSystem, setExpressionSystem] = useState("Homo sapiens");
  const [expressionList, setExpressionList] = useState([]);
  const { PDDispatch } = usePrimerDriverContext();
  const navigate = useNavigate();
  const mode = "PRO";

  useEffect(() => {
    api.data
      .expressionSystem()
      .then(res => setExpressionList(res.data.data))
      .catch(err => console.log(err.message));

    PDDispatch({
      type: "updateMode",
      payload: mode,
    });
  }, [PDDispatch]);

  function validateSequence(sequence) {
    return sequence
      .toUpperCase()
      .split("")
      .filter(char => AMINO_ACIDS.includes(char))
      .join("");
  }

  function handleChangeSequence(e) {
    let { value } = e.target;
    let filteredSequence = validateSequence(value);
    setSequence(filteredSequence);
    setCursorPosition(e.target.selectionStart + 1);
    setSequenceLength(value.length);
  }

  function handleChangeTargetReplacement(value, type) {
    let filteredSequence = validateSequence(value);
    if (type === "target") setTarget(filteredSequence);
    else if (type === "replacement") setReplacement(filteredSequence);
  }

  function handleSubmit(e) {
    e.preventDefault();
    setLoading(true);
    let data;
    api.data
      .primerDriver(formData)
      .then(res => {
        data = res.data;
      })
      .catch(err => {
        console.log(err.message);
        data = "Request failed. Please try again later.";
      })
      .finally(() => {
        PDDispatch({
          type: "updateLoadedResults",
          payload: true,
        });
        PDDispatch({
          type: "updateResults",
          payload: {
            data,
            loaded: true,
          },
        });
        navigate("/results");
      });
  }

  function validateForm() {
    let validSequence = sequenceLength > 0;
    let validMutation = mutationType !== "";
    let validMutationCode =
      mutationType === "sub"
        ? target.length > 0 && replacement.length > 0
        : target.length > 0 || replacement.length > 0;
    setIsValid(validSequence && validMutation && validMutationCode);
    setFormData({
      mode,
      sequence,
      target,
      position,
      replacement,
      mutation_type: mutationType,
      advSettings: {
        Tm_range_min: tmRange[0],
        Tm_range_max: tmRange[1],
        gc_range_min: gcRange[0],
        gc_range_max: gcRange[1],
        length_min: lengthRange[0],
        length_max: lengthRange[1],
        flank5_range_min: flank5Range[0],
        flank5_range_max: flank5Range[1],
        flank3_range_min: flank3Range[0],
        flank3_range_max: flank3Range[1],
        forward_overlap3: forwardOverlap3,
        forward_overlap5: forwardOverlap5,
        terminate_gc: terminateGc,
        center_mutation: centerMutation,
        primer_mode: primerMode,
        expression_system: expressionSystem,
      },
    });
  }

  if (loading) return <LoadingScreen />;
  else
    return (
      <Form
        handleValidate={validateForm}
        handleReset={props.handleReset}
        handleSubmit={handleSubmit}
        isValid={isValid}
        title="Protein-based Primer Design"
      >
        <ProteinSequenceInput
          cursorPosition={cursorPosition}
          handleChangeSequence={handleChangeSequence}
          sequence={sequence}
          sequenceLength={sequenceLength}
        />
        <MutationType handleChange={setMutationType} mutation_type={mutationType} />
        <MutationSelector
          handleChangeTarget={value => handleChangeTargetReplacement(value, "target")}
          handleChangeReplacement={value => handleChangeTargetReplacement(value, "replacement")}
          handleChangePosition={setPosition}
          mutation_type={mutationType}
          position={position}
          replacement={replacement}
          target={target}
        />
        <AdvancedSettings
          tmRange={tmRange}
          setTmRange={setTmRange}
          gcRange={gcRange}
          setGcRange={setGcRange}
          lengthRange={lengthRange}
          setLengthRange={setLengthRange}
          flank3Range={flank3Range}
          setFlank3Range={setFlank3Range}
          flank5Range={flank5Range}
          setFlank5Range={setFlank5Range}
          forwardOverlap3={forwardOverlap3}
          setForwardOverlap3={setForwardOverlap3}
          forwardOverlap5={forwardOverlap5}
          setForwardOverlap5={setForwardOverlap5}
          expressionSystem={expressionSystem}
          setExpressionSystem={setExpressionSystem}
          primerMode={primerMode}
          setPrimerMode={setPrimerMode}
          terminateGc={terminateGc}
          setTerminateGc={setTerminateGc}
          centerMutation={centerMutation}
          setCenterMutation={setCenterMutation}
          expressionList={expressionList}
        />
      </Form>
    );
}

ProteinView.propTypes = {
  handleReset: PropTypes.func.isRequired,
};

export default ProteinView;
