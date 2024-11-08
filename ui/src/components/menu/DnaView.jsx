import PropTypes from "prop-types";
import { useEffect, useState } from "react";
import { useNavigate } from "react-router-dom";
import api from "../../api/index.js";
import { usePrimerDriverContext } from "../../contexts/PrimerDriverContext.jsx";
import {
  AdvancedSettings,
  DnaSequenceInput,
  Form,
  MutationSelector,
  MutationType,
} from "../form/index.js";
import LoadingScreen from "../shared/LoadingScreen.jsx";

const BASES = ["A", "C", "G", "T"];

function DnaView(props) {
  const [formData, setFormData] = useState({});
  const [cursorPosition, setCursorPosition] = useState(1);
  const [sequenceLength, setSequenceLength] = useState(0);
  const [isValid, setIsValid] = useState(false);
  const [loading, setLoading] = useState(false);
  const [sequence, setSequence] = useState("");
  const [mutationType, setMutationType] = useState("");
  const [target, setTarget] = useState("");
  const [replacement, setReplacement] = useState("");
  const [position, setPosition] = useState(0);
  const [tmRange, setTmRange] = useState([75, 85]);
  const [gcRange, setGcRange] = useState([40, 60]);
  const [lengthRange, setLengthRange] = useState([25, 45]);
  const [flank3Range, setFlank3Range] = useState([11, 21]);
  const [flank5Range, setFlank5Range] = useState([11, 21]);
  const [forwardOverlap3, setForwardOverlap3] = useState(9);
  const [forwardOverlap5, setForwardOverlap5] = useState(9);
  const [terminateGc, setTerminateGc] = useState(true);
  const [centerMutation, setCenterMutation] = useState(true);
  const [primerMode, setPrimerMode] = useState("complementary");
  const [expressionSystem, setExpressionSystem] = useState("Homo sapiens");
  const { PDDispatch } = usePrimerDriverContext();
  const navigate = useNavigate();
  const mode = "DNA";

  useEffect(() => {
    PDDispatch({
      type: "updateMode",
      payload: mode,
    });
  }, [PDDispatch]);

  function validateSequence(sequence) {
    return sequence
      .toUpperCase()
      .split("")
      .filter(char => BASES.includes(char))
      .join("");
  }

  function handleChangeSequence(e) {
    const { value } = e.target;
    const filteredSequence = validateSequence(value);
    setSequence(filteredSequence);
    setCursorPosition(e.target.selectionStart + 1);
    setSequenceLength(value.length);
  }

  function handleChangeTargetReplacement(value, type) {
    const filteredSequence = validateSequence(value);
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
        console.error(err.message);
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
    const validSequence = sequenceLength > 0;
    const validMutation = mutationType !== "";
    const validMutationCode =
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
      settings: {
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
  return (
    <Form
      handleValidate={validateForm}
      handleSubmit={handleSubmit}
      handleReset={props.handleReset}
      isValid={isValid}
      title="DNA-based Primer Design"
    >
      <DnaSequenceInput
        cursorPosition={cursorPosition}
        handleChange={handleChangeSequence}
        sequence={sequence}
        sequenceLength={sequenceLength}
      />
      <MutationType handleChange={setMutationType} mutation_type={mutationType} />
      <MutationSelector
        handleChangeTarget={value => handleChangeTargetReplacement(value, "target")}
        handleChangeReplacement={value =>
          handleChangeTargetReplacement(value, "replacement")
        }
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
      />
    </Form>
  );
}

DnaView.propTypes = {
  handleReset: PropTypes.func,
};

export default DnaView;
